get_clmDesign <- function(fullmf, formulas, contrasts) {
### Compute design matrices for a clm object.
### clm-internal+external

    ## Check specification of contrasts:
    checkContrasts(terms=attr(fullmf, "terms"), contrasts=contrasts)

    ## Extract X (design matrix for location effects) + terms, offset:
    res <- get_clmDM(fullmf, formulas, contrasts, type="formula")
    res$off <- res$offset
    res$offset <- NULL

    ## Extract weights:
    res$wts <- getWeights(fullmf)

    ## Extract model response:
    res <- c(get_clmY(fullmf, res$wts), res)

    ## Extract S (design matrix for the scale effects):
    if(!is.null(formulas$scale)) {
        ans <- get_clmDM(fullmf, formulas, contrasts, type="scale")
        res$S <- ans$X
        res$S.terms <- ans$terms
        res$S.off <- ans$offset
        if(attr(ans$terms, "response") != 0)
            stop("response not allowed in 'scale'", call.=FALSE)
    }

    ## Extract NOM (design matrix for the nominal effects):
    if(!is.null(formulas$nominal)) {
        ans <- get_clmDM(fullmf, formulas, contrasts, type="nominal")
        res$NOM <- ans$X
        res$nom.terms <- ans$terms
        if(attr(ans$terms, "response") != 0)
            stop("response not allowed in 'nominal'", call.=FALSE)
        if(!is.null(attr(ans$terms, "offset")))
            stop("offset not allowed in 'nominal'", call.=FALSE)
    }

    ## Return results (list of design matrices etc.):
    res
### NOTE: X, S and NOM are with dimnames and intercepts are
### guaranteed. They may be column rank deficient.
}

get_clmDM <-
    function(fullmf, formulas, contrasts,
             type=c("formula", "scale", "nominal"),
             check.intercept=TRUE, get.offset=TRUE)
{
    ## Get DM (=Design Matrix):
    type <- match.arg(type)
    ## if(type=="scale") browser()
    terms <- terms(formulas[[type]], data=fullmf)
    X <- model.matrix(formulas[[type]], data=fullmf,
                      contrasts.arg=getContrasts(terms, contrasts))
    ## Test for intercept in X(?):
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if(check.intercept && Xint <= 0) {
        X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
        warning(gettext("an intercept is needed and assumed in '%s'", type),
                call.=FALSE)
    } ## Intercept in X is guaranteed.
    res <- list(X=X, terms=terms)
    if(get.offset)
        res$offset <- getOffset(fullmf, terms)
    res
}

get_clmFullmf <- function(mc, fullForm, form.envir, contrasts)
### clm-internal
### Extract the full model.frame
### mc - matched call containing: data, subset, weights, na.action
{
    ## Extract the full model.frame(mf):
    m <- match(c("data", "subset", "weights", "na.action"),
               names(mc), 0)
    mfcall <- mc[c(1, m)]
    mfcall$formula <- fullForm
    mfcall$drop.unused.levels <- TRUE
    mfcall[[1]] <- as.name("model.frame")
    if(is.null(mfcall$data)) mfcall$data <- form.envir
    fullmf <- eval(mfcall, envir=parent.frame(2))
    ## Return:
    fullmf
}

get_clmY <- function(fullmf, wts) {
    y <- model.response(fullmf, "any") ## any storage mode
    if(!is.factor(y)) stop("response needs to be a factor", call.=FALSE)
    ## ylevels are the levels of y with positive weights
    ylevels <- levels(droplevels(y[wts > 0]))
    ## check that y has at least two levels:
    if(length(ylevels) == 1L)
        stop(gettextf("response has only 1 level ('%s'); expecting at least two levels",
                      ylevels), call.=FALSE)
    if(!length(ylevels))
        stop("response should be a factor with at least two levels")
    ## return:
    list(y=y, ylevels=ylevels)
}

get_clmFormulas <- function(mc)
### clm-internal
### Extracts and returns a list of formulas needed for further
### processing.
### mc: matched call
{
    ## Collect all variables in a full formula:
    ## evaluate the formulae in the enviroment in which clm was called
    ## (parent.frame(2)) to get them evaluated properly:
    forms <- list(eval.parent(mc$formula, 2))
    if(!is.null(mc$scale)) forms$scale <- eval.parent(mc$scale, 2)
    if(!is.null(mc$nominal)) forms$nominal <- eval.parent(mc$nominal, 2)
    ## get the environment of the formula. If this does not have an
    ## enviroment (it could be a character), then use the parent frame.
    form.envir <-
        if(!is.null(env <- environment(forms[[1]]))) env
        else parent.frame(2)
    ## ensure formula, scale and nominal are formulas:
    ## forms <- lapply(forms, function(x) {
    ##   try(formula(deparse(x), env = form.envir), silent=TRUE) })
    for(i in 1:length(forms)) {
        forms[[i]] <- try(formula(deparse(forms[[i]]),
                                  env = form.envir), silent=TRUE)
    }
    if(any(sapply(forms, function(f) class(f) == "try-error")))
        stop("unable to interpret 'formula', 'scale' or 'nominal'")
    ## collect all variables in a full formula:
    forms$fullForm <- do.call("getFullForm", forms)
    names(forms)[1] <- "formula"
    ## set environment of 'fullForm' to the environment of 'formula':
    forms$form.envir <- environment(forms$fullForm) <- form.envir
    ## return:
    forms
}

get_clmRho <-
    function(fullmf, formulas, contrasts, link, threshold,
             parent=parent.frame(), start=NULL, ...)
### .default method(?)
{
    ## Get design matrices etc:
    design <- get_clmDesign(fullmf=fullmf,
                            formulas=formulas,
                            contrasts=contrasts)
    ## Get threshold information:
    design$ths <- makeThresholds(design$ylevels, threshold)
    ## Drop columns for aliased coefs:
    design <- drop.cols(design, drop.scale=FALSE, silent=TRUE)
    ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
    rho <- with(design, {
        clm.newRho(parent.frame(), y=y, X=X, NOM=design$NOM, S=design$S,
                   weights=wts, offset=off, S.offset=design$S.off,
                   tJac=ths$tJac)
    })
    ## Set and check starting values for the parameters:
    start <- set.start(rho, start=start, get.start=is.null(start),
                       threshold=threshold, link=link, frames=design)
    rho$par <- as.vector(start) ## remove attributes
    ## Set pfun, dfun and gfun in rho:
    setLinks(rho, link)
    ## Return:
    rho
}

get_clmRho.clm <-
    function(object, parent=parent.frame(), ...) {
### Safely generate the model environment from a model object.
    o <- object
    get_clmRho(fullmf=model.frame(o), formulas=o$formulas,
               contrasts=o$contrasts, start=c(o$start), link=o$link,
               threshold=o$threshold, parent=parent, ...)
}

## get_mfcall <- function(mc, envir=parent.frame(2)) {
##     m <- match(c("data", "subset", "weights", "na.action"),
##                names(mc), 0)
##     mf <- mc[c(1, m)]
##     ## mf$formula <- fullForm
##     mf$drop.unused.levels <- TRUE
##     mf[[1]] <- as.name("model.frame")
##     ## if(is.null(mf$data)) mf$data <- form.envir
##     list(mfcall=mf, envir=parent.frame(2))
## }

## get_clmDesign <- function(fullmf, formulas, contrasts) {
## ### Compute design matrices for a clm object.
## ### clm-internal+external
##
##     ## Extract X (design matrix for location effects) + terms, offset:
##     res <- get_clmX(fullmf, formulas, contrasts)
##
##     ## Extract weights:
##     res$wts <- getWeights(fullmf)
##
##     ## Extract model response:
##     res <- c(get_clmY(fullmf, res$wts), res)
##
##     ## Extract S (design matrix for the scale effects):
##     if(!is.null(formulas$scale))
##         res <- c(res, get_clmS(fullmf, formulas, contrasts))
##
##     ## Extract NOM (design matrix for the nominal effects):
##     if(!is.null(formulas$nominal))
##         res <- c(res, get_clmNOM(fullmf, formulas, contrasts))
##     ## Return results (list of design matrices etc.):
##     res
## ### NOTE: X, S and NOM are with dimnames and intercepts are
## ### guaranteed. They may be column rank deficient.
## }


## get_clmX <- function(fullmf, formulas, contrasts) {
##     ## mfcall$formula <- formulas$formula
##     ## mfcall$na.action <- "na.pass" ## filter NAs by hand below
##     ## X.mf <- eval(mfcall, envir=mfenvir)
##     ## X.terms <- attr(X.mf, "terms")
##     ## X <- model.matrix(X.terms, data=fullmf,
##     ##                 contrasts.arg=getContrasts(X.terms, contrasts))
## ### NOTE: The following construct is not enough, since calls like
## ### clm(ordered(rating) ~ temp + contact, data=wine) will fail. The
## ### reason is that 'ordered(factor)' will not be evaluated correctly.
## ### X.mf <- model.frame(formulas$formula, data=mf)
## ### mf$formula <- forms[[1]]
## ### X.mf <- eval(mf, envir = parent.frame(2))
## ### X.terms <- attr(X.mf, "terms")
## ### X <- model.matrix(X.terms, data=mf,
## ###                   contrasts.arg=getContrasts(X.terms, contrasts))
## ###
## ### Can we somehow extract the right colums/terms from mf anyway?
## ###
##     terms <- terms(formulas$formula, data=fullmf)
## ### FIXME: test if fm$terms contains contact (which is should not) if
## ### fm <- clm(rating ~ terms, scale=~contact, data=wine)
##     X <- model.matrix(formulas$formula, data=fullmf,
##                       contrasts.arg=getContrasts(terms, contrasts))
##     n <- nrow(X)
##     ## Test for intercept in X:
##     Xint <- match("(Intercept)", colnames(X), nomatch = 0)
##     if(Xint <= 0) {
##         X <- cbind("(Intercept)" = rep(1, n), X)
##         warning("an intercept is needed and assumed in 'formula'",
##                 call.=FALSE)
##     } ## intercept in X is guaranteed.
##     off <- getOffset(fullmf, terms)
##     ## Filter NAs if any:
##     ## if(!is.null(naa <- na.action(fullmf)))
##     ##     off <- off[-naa]
##     ## return:
##     list(X=X, terms=terms, off=off)
## }
##
## get_clmS <- function(mfcall, mfenvir, fullmf, formulas, contrasts) {
##     mfcall$formula <- formulas$scale
##     mfcall$na.action <- "na.pass" ## filter NAs by hand below
##     mf <- eval(mfcall, envir=mfenvir)
##     terms <- attr(mf, "terms")
##     X <- model.matrix(terms, data=fullmf,
##                       contrasts.arg=getContrasts(terms, contrasts))
##     ## S.mf <- model.frame(formulas$scale, data=fullmf)
##     if(!is.null(model.response(mf)))
##         stop("response not allowed in 'scale'", call.=FALSE)
##     ## S.terms <- attr(mf, "terms")
##     ## S <- model.matrix(S.terms, data=fullmf,
##     ##                   contrasts.arg=getContrasts(S.terms, contrasts))
##     ## Test for intercept in S:
##     Xint <- match("(Intercept)", colnames(X), nomatch = 0)
##     if(Xint <= 0) {
##         X <- cbind("(Intercept)" = rep(1, n), X)
##         warning("an intercept is needed and assumed in 'scale'",
##                 call.=FALSE)
##     } ## intercept in S is guaranteed.
##     off <- getOffset(mf, terms)
##     ## Filter NAs if any:
##     if(!is.null(naa <- na.action(fullmf)))
##         off <- off[-naa]
##     ## Return:
##     list(S=X, S.terms=terms, S.off=off)
## }
##
## get_clmNOM <- function(mfcall, mfenvir, fullmf, formulas, contrasts) {
##     ## ans <- get_clmDM(mfcall, mfenvir, fullmf, formulas, contrasts,
##     ##                  type="nominal", get.offset=FALSE)
##     ## if(!is.null(model.response(ans$mf)))
##     ##     stop("response not allowed in 'nominal'", call.=FALSE)
##     ## if(!is.null(model.offset(ans$mf)))
##     ##     stop("offset not allowed in 'nominal'", call.=FALSE)
##     ## return(list(NOM=ans$X, nom.terms=terms))
##     mfcall$formula <- formulas$nominal
##     mfcall$na.action <- "na.pass" ## filter NAs by hand below
##     mf <- eval(mfcall, envir=mfenvir)
##     terms <- attr(mf, "terms")
##     X <- model.matrix(terms, data=fullmf,
##                       contrasts.arg=getContrasts(terms, contrasts))
##     ## nom.mf <- model.frame(formulas$nominal, data=fullmf)
##     if(!is.null(model.response(mf)))
##         stop("response not allowed in 'nominal'", call.=FALSE)
##     if(!is.null(model.offset(mf)))
##         stop("offset not allowed in 'nominal'", call.=FALSE)
##     ## terms <- attr(nom.mf, "terms")
##     ## NOM <- model.matrix(terms, data=fullmf,
##     ##                     contrasts.arg=getContrasts(terms, contrasts))
##     Xint <- match("(Intercept)", colnames(X), nomatch = 0)
##     if(Xint <= 0) {
##         X <- cbind("(Intercept)" = rep(1, n), X)
##         warning("an intercept is needed and assumed in 'nominal'",
##                 call.=FALSE)
##     } ## intercept in NOM is guarantied.
##     list(NOM=X, nom.terms=terms)
## }
