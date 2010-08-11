drop.coef <- function(X, silent = FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0 
{
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
  if(qr.X$rank == ncol(X))
    return(X) ## return X if X has full column rank
  if(!silent) ## message the no. dropped columns:
    message(gettextf("design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  ## did we succeed? stop-if-not:
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX) 
}

