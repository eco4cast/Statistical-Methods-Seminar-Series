
## copied from glmmTMB:
reOnly <- function(f,response=FALSE,bracket=TRUE) {
    ff <- f
    if (bracket)
        ff <- lapply(lme4::findbars(ff),makeOp,quote(`(`)) ## bracket-protect terms
    ff <- Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),ff)
    if (response && length(f)==3) {
        form <- makeOp(f[[2]],ff,quote(`~`))
    } else {
        form <- makeOp(ff,quote(`~`))
    }
    return(form)
}

sumTerms <- function(termList) {
    Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),termList)
}

## combine unary or binary operator + arguments (sugar for 'substitute')
## FIXME: would be nice to have multiple dispatch, so
## (arg,op) gave unary, (arg,arg,op) gave binary operator
makeOp <- function(x,y,op=NULL) {
    if (is.null(op)) {  ## unary
        substitute(OP(X),list(X=x,OP=y))
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

g4fit <- function(formula,data,na.action,...) {
    ## remove all NA-containing rows
    mf <- data[all.vars(formula)]
    na.removed <- attr(na.action(mf),"na.action")
    data <- data[complete.cases(mf),]
    ## now fit gamm4
    g <- gamm4(formula=nobars(formula),
               random=as.formula(reOnly(formula)),
               data=data,...)
    ## now restore information
    attr(g$mer@frame,"na.action") <- na.removed
    class(g) <- c("gamm4","list")
    return(g)
}

logLik.gamm4 <- function(object, ...) {
    logLik(object$mer, ...)
}

tidy.gamm4 <- function(x,...) {
    r <- tidy(x$mer,...)
    r$term <- gsub("^X","",r$term)
    return(r)
}

## ugh. result of trial-and-error with object_size(), neutralizing
##  environments where we can find them.  Wonder if extractors still
##  work?
strip_gamm4_env <- function(x) {
    gam_env_obj <- c("pred.formula","formula","pterms","terms","model")
    for (i in gam_env_obj) {
        environment(x$gam[[i]]) <- NULL
    }
    x$gam$model <- c(x$gam$model)
    mer_env_obj <- c("pp","frame")
    for (i in mer_env_obj) {
        environment(slot(x$mer,i)) <- NULL
    }
    ## should strip environment from formula; don't trash formula ...
    formula <- attr(x$mer@frame,"formula")
    environment(formula) <- NULL
    terms <- attr(x$mer@frame,"terms")
    environment(terms) <- NULL
    x$mer@frame <- as.data.frame(c(x$mer@frame))
    attr(x$mer@frame,"formula") <- formula
    attr(x$mer@frame,"terms") <- terms
    return(x)
}

formula.gamm4 <- function(x,fixed.only=FALSE) {
    if (fixed.only) return(formula(x$gam))
}

## https://stats.stackexchange.com/questions/65643/using-a-gamm4-model-to-predict-estimates-in-new-data
predict.gamm4 <- function(x,re.form=NULL,newdata=NULL,...) {
    if (is.null(newdata) && is.null(re.form)) {
        ## the easy case: just extracts fitted values
        return(predict(x$mer,re.form=re.form, ...))
    }
    newdata <- data.frame(newdata,Xr=NA)  ## fake column for smooth spline term
    ## if re.form=NULL (all REs) specified, need to drop 1|Xr term from REs
    fb <- lme4::findbars(formula(x$mer))
    XrTerm <- sapply(fb,function(x) identical(x,quote(1|Xr)))
    if (is.null(re.form)) {
        re.form <- sumTerms(fb[-XrTerm])
        re.form <- as.formula(makeOp(re.form,quote(`~`)))
    }
    ## predict.gam has fixed effects + spatial smooths
    fixed <- predict(x$gam,newdata=newdata,...)
    ## predict.random gets *only* non-spatial RE terms (0 if none specified)
    random <- predict(x$mer,newdata=newdata,random.only=TRUE,re.form=re.form,...)
    return(fixed+random)
}

terms.gamm4 <- function(x) {
    terms(x$gam)
}

formula.gamm4 <- function(x,type=c("gam","mer"),drop.smooth=TRUE,fixed.only=FALSE,...) {
    type <- match.arg(type)
    ff <- formula(x[[type]])
    if (drop.smooth) {
        ff <- glmmTMB::drop.special(ff,quote(s))
    }
    return(ff)
}


vcov.gamm4 <- function(x,type=c("mer","gam")) {
    type <- match.arg(type)
    vcov(x[[type]])
}

coef.gamm4 <- function(x,type=c("mer","gam")) {
    type <- match.arg(type)
    coef(x[[type]])
}

fixef.gamm4 <- function(x) {
    ff <- fixef(x$mer)
    names(ff) <- gsub("^X","",names(ff))
    return(ff)
}

ranef.gamm4 <- function(x) {
    ranef(x$mer)
}

