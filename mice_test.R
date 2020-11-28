#---- simulate data
library(mice)
library(miceadds)


dummy <- function(){
  set.seed(987)
  N <- 1000         # number of persons
  var_err <- .4     # error variance
  dat <- data.frame( x1=stats::rnorm(N), x2=stats::rnorm(N) )
  dat$theta <- .3 * dat$x1 - .5*dat$x2 + stats::rnorm(N)
  dat$y <- dat$theta + stats::rnorm( N, sd=sqrt(var_err) )
  browser()
  #-- linear regression for measurement-error-free data
  mod0a <- stats::lm( theta ~ x1 + x2, data=dat )
  summary(mod0a)
  #-- linear regression for data with measurement error
  mod0b <- stats::lm( y ~ x1 + x2, data=dat )
  summary(mod0b)
  
  #-- process data for imputation
  
  dat1 <- dat
  dat1$theta <- NA
  scale.values <- list( "theta"=list( "M"=dat$y, "SE"=rep(sqrt(var_err),N )))
  dat1$y <- NULL
  
  cn <- colnames(dat1)
  V <- length(cn)
  method <- rep("", length(cn) )
  names(method) <- cn
  method["theta"] <- "plausible.values"
  
  #-- imputation in mice
  imp <- mice::mice( dat1, maxit=1, m=5, allow.na=TRUE, method=method,
                     scale.values=scale.values )
  summary(imp)
  
  #-- inspect first dataset
  summary( mice::complete(imp, action=1) )
  
  #-- linear regression based on imputed datasets
  mod1 <- with(imp, stats::lm( theta ~ x1 + x2 ) )
  summary( mice::pool(mod1) )
  
  ## End(Not run)
  
}
