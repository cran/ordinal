## This file contains:
## Functions (getThetamat) to compute a table of threshold
## coefficients from model fits (clm()s) with nominal effects.

getThetamat <-
  function(terms, alpha, assign, contrasts, xlevels, tJac, dataClasses)
### Compute matrix of thresholds for all combinations of levels of
### factors in the nominal formula.
###
### Input:
### terms: nominal terms object
### alpha: vector of threshold parameters
### assign: attr(NOM, "assign"), where NOM is the design matrix for
###   the nominal effects
### contrasts: list of contrasts for the nominal effects
### tJac: threshold Jacobian with appropriate dimnames.
###
### Output:
### Theta: data.frame of thresholds
### mf.basic: if nrow(Theta) > 1 a data.frame with factors in columns
###   and all combinations of the factor levels in rows.
{
    ## Make matrix of thresholds; Theta:
    Theta <- matrix(alpha, ncol=ncol(tJac), byrow=TRUE)
    ## matrix with variables-by-terms:
    factor.table <- attr(terms, "factors")
    all.varnm <- rownames(factor.table)
### NOTE: need to index with all.varnm not to include (weights) and
### possibly other stuff.
    var.classes <- dataClasses[all.varnm]
  numeric.var <- which(var.classes != "factor")
### NOTE: for "dataClasses" see help(.MFclass). logical variables are
### treated like numeric variables.
  ## terms associated with numerical variables:
    numeric.terms <- factor.terms <- numeric(0)
    if(length(factor.table)) {
        numeric.terms <-
            which(colSums(factor.table[numeric.var, , drop=FALSE]) > 0)
        ## terms only involving factor variables:
        factor.terms <-
            which(colSums(factor.table[numeric.var, , drop=FALSE]) == 0)
    }
    ## remove rows in Theta corresponding to numerical variables:
    if(length(numeric.terms)) {
### NOTE: ncol(NOM) == length(asgn) == nrow(Theta)
### length(attr(terms, "term.labels")) == ncol(factor.table)
### NOTE: length(var.classes) == nrow(factor.table)
        numeric.rows <- which(assign %in% numeric.terms)
        Theta <- Theta[-numeric.rows, , drop=FALSE]
        if(length(factor.terms))
            terms <- drop.terms(terms, dropx=numeric.terms,
                                keep.response=FALSE)
    }
    ## if some nominal effects are factors:
    if(length(factor.terms)) {
        ## get xlevels for factors, not ordered (factors)
        factor.var <- which(var.classes == "factor")
        factor.varnm <- names(var.classes)[factor.var]
        xlev <- xlevels[factor.varnm]
        ## minimal complete model frame:
        mf.basic <- do.call(expand.grid, xlev)
        ## minimal complete design matrix:
        X <- model.matrix(terms, data=mf.basic,
                          contrasts=contrasts)
### NOTE: There are no contrasts for numerical variables.
### FIXME: remove contrasts for 'ordered' variables.
        ## from threshold parameters to thresholds:
        Theta <- apply(Theta, 2, function(th) X %*% th)
    }
    ## adjust each row in Theta for threshold functions:
    tmp <- lapply(seq_len(nrow(Theta)), function(i)
                  c(tJac %*% Theta[i, ]))
    Theta <- do.call(rbind, tmp)
### NOTE: apply returns a vector and not a matrix when ncol(Theta) ==
### 1, so we need to avoid it here.
    ## Theta <- t(apply(Theta, 1, function(th) tJac %*% th))
    colnames(Theta) <- rownames(tJac)
    res <- list(Theta = as.data.frame(Theta))
    ## add factor information if any:
    if(nrow(Theta) > 1)  res$mf.basic <- mf.basic
    ## return:
    res
}
