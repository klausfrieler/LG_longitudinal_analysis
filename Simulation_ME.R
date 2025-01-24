### ME simulations
##Psychoacoustic Scenario
# => In the below all variables are linearly d and z-transformed;
# It would be more realistic if y was expressed in Hz, with a skewed distributions ranging from 330 (reference frequency) and 400
#set.seed(10081971)
require(univOutl)
require(lavaan)
require(TSI) # install from github: install.packages("remotes") /n remotes::install_github("mmansolf/TSI")
library(tidyverse)

messagef <- function(...) message(sprintf(...))
scale_ <- function(x) as.numeric(scale(x))

generate_data <- function(b0, b1, b2, error_type, me = 1, n = 50, n_batch = 50, normal = F){
  n <- n * n_batch
  if(!normal){
  xt <- scale_(runif(
    n = n,
    min = -1,
    max = 1))
  } else{
    xt <- rnorm(n)
  }
  
  if(!normal){
    zt <- scale_(18 + rpois(n = n, lambda = 5))
  } 
  else{
    zt <- rnorm(n, 0, 5)
  }
  
  yt <- b0 + b1 * xt + b2 * zt
  ye <- rnorm(n = n,
              mean = 0,
              sd = me * sd(yt))
  
  if (error_type == "classical") {
    y <- yt + ye
    #print(colnames(df))
  }
  else if (error_type == "systematic") {
    y <- me + (1 + me) * yt + ye
    #browser()
  }
  else if (error_type == "heteroscedastic") {
    ye <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * yt^2))
    y <- yt + ye
  }
  else if (error_type == "differential") {
    ye <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * xt^2))
    y <- yt + ye
    #y <- me + (1 + me) * xt + ye + yt
  }
  else{
    stop(sprinf("Unknow error type: %s", error_type))
  }
  tibble(y = y, yt = yt, ye = ye, xt = xt, zt = zt, batch = rep(1:(n/n_batch), n_batch))
  
}

get_coefs <- function(df, method, measurement_error){
  if (method == "no_correction") {
    coefs <- coef(lm(y ~ xt + zt, data = df))
  }
  else if (method == "outlier_exclusion") {
    ol <- suppressMessages(
      boxB(x = df$y,
           k = 1.5,
           method = 'resistant')$outliers)
    if (length(ol) >= 1) {
      coefs <- coef(lm(y ~ xt + zt, data = df[-ol, ]))
    }
    else{
      coefs <- coef(lm(y ~ xt + zt, data = df))
    }
  }
  
  # if (methods[j] == "weighting") {
  #   inv_err <- rep(mean(1 / df$ye ^ 2), nrow(df))
  #   coefs <- coef(lm(y ~ xt + zt, weights = inv_err, data = df))
  #   eval_w <- rbind(eval_w, eval(
  #     coefs = coefs,
  #     b0 = b0,
  #     b1 = b1,
  #     b2 = b2
  #   ))
  # }
  
  else if (method == "LV") {
    r.ly <- 1 / (1 + measurement_error ^ 2) #compute reliability of latent variable (i.e. commn variance) from measurement error
    r.y <- 1 - r.ly    #compute unique variance from measurement error
    m <- '
        ly =~ 1 * y
        ly ~~ r.ly * ly
        y ~~ r.y * y
        ly ~ xt + zt
        '
    fit <- sem(
      m,
      data = df,
      estimator = "MLR",
      meanstructure = T
    )
    coefs <- parameterEstimates(fit)[c(9, 4, 5), "est"]
  }
  else if (method == "MI") {
    #browser()
    mice.dfs = TSI(
      df %>% select(-c(yt, ye)),
      os_names = c("y"),
      #se_names = c("ye"),
      metrics = "z",
      score_types = "CTT",
      reliability = c(1 / (1 + measurement_error ^ 2)),
      separated = F,
      mice_args = c(
        m = 5,
        maxit = 5,
        printFlag = F
      )
    )
    #browser()
    coefs <- pool(with(mice.dfs, lm(y ~ xt + zt))) %>% pluck("pooled") %>% pluck("estimate")
  } 
  else{
    stop(sprintf("Unknown method: %s", method))
  }
  coefs
}

eval <- function(coefs, true_coefs) {
  #browser()
  names(coefs) <- names(true_coefs)
  
  # rbias_b0 <- abs(b0 - coefs[1]) / abs(b0)
  # rbias_b1 <- abs(b1 - coefs[2]) / abs(b1)
  # rbias_b2 <- abs(b2 - coefs[3]) / abs(b2)
  #bias <-  c(rbias_b0, rbias_b1, rbias_b2)
  abs_error <- abs(coefs- true_coefs)/abs(true_coefs)
  rel_error <- abs_error/abs(true_coefs)
  t(rbind(true_coefs, coefs, abs_error, rel_error)) %>% as.data.frame() %>% rownames_to_column("name") %>% as_tibble()
}

gen_smp_sys_error_1 <- function(b0 = 1,
                                b1 = -0.5,
                                b2 = 0.25,
                                error_type = c("classical", "systematic", "heteroscedastic", "differential"),
                                measurement_error = c(1),
                                n = 50,
                                n_batch = 50,
                                normal = TRUE,
                                methods = c("no_correction", 
                                            "outlier_exclusion", 
                                            #"Weighting", 
                                            "LV",
                                            "MI")) {
  #the amount of measure error (me) is a parameter, possible values are: 0, 0.33, 0.65, 1, 1.52
  #sample size (n) is a parameter, possible values are 50, 500, 5000
  #include parameter for level of ME info: individual-level vs. variable-level?
  #include parameter for true ME info vs. aggregated/approximate ME info?
  #error_type  <- match.arg(error_type) 
  map_dfr(measurement_error, function(me){
    map_dfr(error_type, function(et){
      messagef("Simulating data with '%s' error (me = %.2f)...", et, me)
      simulated_data <- generate_data(b0, b1, b2, error_type = et,
                                      me = me,
                                      n = n, n_batch = n_batch, 
                                      normal = normal)
      
      map_dfr(1:n_batch, function(i){
        #simulation
        #browser()
        #analysis methods and evaluation, i.e. comparison of coefficients to ground truth from simulation
        df <- simulated_data[simulated_data$batch == i, ] %>% select(-batch)
        map_dfr(methods, function(method){
          #browser()
          coefs <- get_coefs(df, method, me)
          eval(
            coefs, c(
              "b0" = b0,
              "b1" = b1,
              "b2" = b2
            )) %>%
            mutate(method = !!method, 
                   batch = i, 
                   error_type = et)
        })
      })
    }) %>% 
      mutate(n = n, 
             n_batch = n_batch, 
             measurement_error = me)
  })
}

