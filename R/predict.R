##########################################################################
## predict method for mnlogit objects                                    #
## Contributed by: Florian Oswald, University College London             #
## Homepage: http://floswald.github.io					 #
##########################################################################
predict.mnlogit <- function(object, newdata=NULL, probability=TRUE,
                            returnData = FALSE, choiceVar=NULL, ...) 
{
    size     <- object$model.size
    # get choice set for colnames
    choiceSet <- unique(index(object)$alt)
    #kk: choiceSet <- unique(newdata$LCagg_model)  #xavi: needs to be generalised!!
    if (is.null(newdata)) {
        # if no new data, use probabilities computed during training model
        if (probability)
	          return(object$probabilities)
        else { 
	          return(apply(object$probabilities, 1, function(x)
			                  object$choices[which(x == max(x, na.rm = TRUE))]))
        }
    } else {
	      # make sure newdata is ordered by choice
        if (is.null(choiceVar)) {
          if (!any(class(newdata) == "mlogit.data"))
            stop("NULL choiceVar requires newdata to be a mlogit.data object")
          if (nrow(newdata) != nrow(attr(newdata, "index")))
            stop("mlogit.data object newdata has incorrect index attribute")
          choiceVar <- "_Alt_Indx_"
          newdata[[choiceVar]] <- attr(newdata, "index")$alt
          #newdata[[choiceVar]] <- index(object)$alt
        }
	      newdata <- newdata[order(newdata[[choiceVar]]), ]

	      # Get name of response column
	      pf <- parseFormula(object$formula)
	      resp.col <- attr(pf, "response")   #xavi: the one with 0 and 1

        # check that all columns from data are present (except response col)
        # this is important when you build Y below.
	      newn <- names(newdata)
	      oldn <- setdiff(names(object$model), resp.col)
      	if (!all(oldn %in% newn))
      	    stop("newdata must have same columns as training data. ")

	      # different model size: N # newdata must have N*K rows
	      # size$K <- length(unique(newdata[[choiceVar]])) #xavi: Modifying K, number of choices of newdata given that could be different from the ones of the model fitted
	                                               #xavi: This might happen with choices given in few cases. They might not be present when split data for calibration and prediction
	      
	      if (nrow(newdata) %% size$K)   #xavi: if number of choices in newdata and in fitting data are different, stop it
	                                     #xavi: avoid this point changing size$K with number of choices in newdata 
	        #xavi: stop("Mismatch between nrows in newdata and number of choices.")
	        print("Mismatch between nrows in newdata and number of choices. # xavi: however, go ahead!!")
    }
    data <- newdata
    size$N <- nrow(data)/size$K       # number of individuals #xavi: in newdata
    if (!(resp.col %in% names(data))) # attach a response column 
        data[[resp.col]] <- rep(1, size$N)

    # Initialize utility matrix: dim(U) = N x K-1
    probMat <- matrix(rep(0, size$N * (size$K-1)), nrow=size$N, ncol=size$K-1)

    formDesignMat <- function(varVec = NULL, includeIntercept = TRUE)
    {
        if (is.null(varVec) && !includeIntercept) return(NULL) 
        fm <- paste(attr(formula, "response"), "~")
        if (!is.null(varVec))
            fm <- paste(fm, paste(varVec, collapse = "+")) #xavi:this is the formula
        if (!includeIntercept) fm <- paste(fm, "-1 ")
        else fm <- paste(fm, "+1 ")
        modMat <- model.matrix(as.formula(fm), data)  #xavi: made with newdata
    }
    # Grab the parsed formula from the fitted mnlogit object 
    formula  <- parseFormula(object$formula)

    X <- formDesignMat(varVec = attr(formula, "indSpVar"), 
                       includeIntercept = attr(formula, "Intercept"))
    X <- if (!is.null(X)) X[1:size$N, , drop=FALSE]   # Matrix of ind sp vars #xavi: we need to keep also, and only!, the variables not removed during calibration because of their linear dependency
    
    Y <- formDesignMat(varVec = attr(formula, "csvChCoeff"), 
                       includeIntercept = FALSE)
    Z <- formDesignMat(varVec = attr(formula, "csvGenCoeff"), 
                       includeIntercept = FALSE)

    # Do the subtraction: Z_ik - Zi0 (for Generic coefficients data)
    ### NOTE: Base choice (with respect to normalization) is fixed here
    ###       Base choice is the FIRST alternative
    if(!is.null(Z)) { 
        for (ch_k in 2:size$K) {
            Z[((ch_k - 1)*size$N + 1):(ch_k*size$N), ] <-
              Z[((ch_k-1)*size$N+1):(ch_k*size$N), , drop=FALSE] - 
	            Z[1:size$N, , drop=FALSE]
        }
    }
    # Drop rows for base alternative
    Z <- Z[(size$N + 1):(size$K*size$N), , drop=FALSE]
    
    #xavi: Remove variables that have been removed during calibration because of linear dependency
    if (!is.null(X)){
      indSpVar2keep <- object$varNamesList$indSpVar
      indSpVar2keep <- indSpVar2keep
      X <- X[, indSpVar2keep]
    }

    if (!is.null(Y)){
      csvChCoeffVar2keep <- object$varNamesList$csvChCoeff
      csvChCoeffVar2keep <- csvChCoeffVar2keep
      Y <- Y[, csvChCoeffVar2keep]
    }
    
    if (!is.null(Z)){
      csvGenCoeffVar2keep <- object$varNamesList$csvGenCoeff
      csvGenCoeffVar2keep <- csvGenCoeffVar2keep
      Z <- Z[, csvGenCoeffVar2keep]
    }
    
    # Grab trained model coeffs from fitted mnlogit object
    coeffVec <- object$coeff

    # First compute the utility matrix (stored in probMat)
    if (size$p) {  #xavi: this is for individual-specific variables
         probMat <- probMat + X %*% matrix(coeffVec[1:((size$K-1) *size$p)],
			                          nrow = size$p, ncol = (size$K-1), byrow=FALSE)
    }
    if (size$f) {
        findYutil<- function(ch_k)
	      {
	          offset <- (size$K - 1)*size$p
	          init <- (ch_k - 1)*size$N + 1
	          fin <- ch_k * size$N
	          Y[init:fin, , drop=FALSE] %*%
	          coeffVec[((ch_k-1)*size$f + 1 + offset):(ch_k*size$f+offset)]
	      }
	      vec <- as.vector(sapply(c(1:size$K), findYutil))
	      # normalize w.r.t. to k0 here - see vignette on utility normalization
	      vec <- vec - vec[1:size$N]
	      probMat <- probMat + matrix(vec[(size$N+1):(size$N*size$K)], 
				                         nrow = size$N, ncol = (size$K-1), byrow=FALSE)
    }
    if (size$d) {
	      probMat <- probMat +
	          matrix(Z %*% coeffVec[(size$nparams - size$d + 1):size$nparams],
		            nrow = size$N, ncol=(size$K-1), byrow=FALSE)
    }

    # Convert utility to probabilities - use logit formula
    probMat <- exp(probMat)                           # exp(utility)

    #xavi: baseProbVec <- 1/(1 + rowSums(probMat))           # P_i0
    baseProbVec <- 1/(1 + rowSums(probMat, na.rm = TRUE))   # P_i0  #xavi: to avoid NA as results when NA present in data frame
    probMat <- probMat * matrix(rep(baseProbVec, size$K-1),
		      nrow = size$N, ncol = size$K-1) # P_ik

    probMat <- cbind(baseProbVec,probMat)

    if (nrow(probMat) == 1)
	      probMat <- as.matrix(probMat)
    colnames(probMat) <- choiceSet

    if (probability) { #xavi:this gives probs
         if (returnData) attr(probMat, "data") <- newdata
	      return(probMat)
    } else { #this gives the best choice
	      choice <- apply(probMat, 1, function(x)
			                  object$choices[which(x == max(x, na.rm = TRUE))])
        if (returnData) attr(choice, "data") <- newdata
	      return(choice)
    }
}
