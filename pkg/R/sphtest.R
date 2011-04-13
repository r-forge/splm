
sphtest <- function (x, ...)
{
    UseMethod("sphtest")
}

sphtest.formula <- function (x, data, w, lag=F, alternative = c("random", "pooled"), ...)
{
 ## performs a Hausman test of a FE model with spatial lag or error
 ## against "alternative" with same spatial specification
    errors <- switch(match.arg(alternative),
			random = {if(lag) "re" else "semre"},
			pooled = {if(lag) "ols" else "sem"})

    if("listw" %in% class(w)) {
      W <- listw2mat(w)
      lw <- w
    } else {
      W <- w
      lw <- mat2listw(w)
    }

    model <- if(lag) "lag" else "err"

    x0 <- update(x, .~.-1)
    
    femod <- spfeml(x0, data, listw = lw,
                    model = model, effects = "spfe",
                    ...)

    remod <-spreml(x, data, w = W, errors=errors,
                   lag=lag, ...)
    
    sphtest(femod, remod, ...)
}

sphtest.splm <- function (x, x2, ...)
{
  ## check if lag or error
    is.lag <- function(mod) {
      if(is.null(mod$spat.coef)) {
        ## model is random or pooled
        is.lag <- !(is.null(mod$arcoef))
      } else {
        ## model is FE
        is.lag <- names(mod$spat.coef)=="rho"
      }
      return(is.lag)
    }
  ## check consistency
    if(is.lag(x)!=is.lag(x2)) stop("Models are heterogeneous")
    
  ## test on coefficients (excluding SAR)      
  ## model order is irrelevant
    coef.wi <- coef(x)[-1]
    coef.re <- coef(x2)[-1]
    vcov.wi <- x$vcov[-1,-1]
    vcov.re <- x2$vcov[-1,-1]
    names.wi <- names(coef.wi)
    names.re <- names(coef.re)
    dbeta <- coef.wi - coef.re
    df <- length(dbeta)
    dvcov <- vcov.re - vcov.wi
    stat <- abs(t(dbeta) %*% solve(dvcov) %*% dbeta)
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    names(stat) <- "chisq"
    parameter <- df
    names(parameter) <- "df"
    data.name <- paste(deparse(x$call$formula))
    alternative <- "one model is inconsistent"
    res <- list(statistic = stat, p.value = pval, parameter = parameter,
        method = "Hausman test for spatial models", data.name = data.name, alternative = alternative)
    class(res) <- "htest"
    return(res)
}
