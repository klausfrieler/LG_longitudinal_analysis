library(tidyverse)
library(log4r)
library(tictoc)
logger <- create.logger()
logfile(logger) <- "me_simu.log"
level(logger) <- "INFO"

version <- "0.4.0"
aws3_dir <- "measurement_error_study"

bar <- paste(rep("-", 60), collapse ="")

#This hack avoids including library(mclust), due to bug in mclust
mclustBIC <- mclust::mclustBIC

source("setup.R")
aws.set_credentials("kf_aws-credentials.txt")

messagef <- function(...){
  message(sprintf(...))
  log4r::info(logger, sprintf(...))
}

deco_messagef <- function(...){
  messagef(bar)
  messagef(...)
  messagef(bar)
}
#source("scripts/IRT_bootstrapper.R")

required_packages <- c("aws.s3", "lubridate", 
                       "brms", "mice", "miceadds", "mclust", 
                       "catR", "Amelia", "broom", "MASS", 
                       "faux", "broom.mixed", "simex", "digest")

prepare_simulations <- function(){
  for(n in required_packages){
    installed <- n %in% installed.packages()
    messagef("Checking for package %s: %s",n, ifelse(installed, "Installed", "Not installed"))
  }
  setup_workspace_me(F)
  messagef("Running test simulation... bear with me.")
  tmp <- test_simulations(master_cross, n_simul = 1)
}

imputed_lm <- function(imputed, model = lm_model, alpha = .05){
  #browser()
  tidied <- 
    map_dfr(1:length(imputed$imputations), function(i){
      lm_wrapper(imputed$imputations[[i]], broom::tidy) %>% mutate(iter = i)
  })
  glanced <-  
    map_dfr(1:length(imputed$imputations), function(i){
      lm_wrapper(imputed$imputations[[i]], broom::glance) %>% mutate(iter = i)
  })
  
  qm <- quick_meld(tidied, by_term = T) %>% 
    mutate(df.residual = glanced$df.residual[1], adj.r.squared = mean(glanced$adj.r.squared)) %>% 
    mutate(statistic = estimate / std.error,
           low_ci = estimate + std.error * qt(alpha/2, df.residual),
           up_ci = estimate + std.error * qt(1-alpha/2, df.residual),
           p.value = 2 * pt(abs(statistic), df.residual, lower.tail = FALSE)) 
  qm
}

quick_meld <- function(imputed_model, by_term = F){
  if(by_term){
    map_dfr(unique(imputed_model$term), function(t){
      tmp <- imputed_model %>% filter(term == t)
      pool <- Amelia::mi.meld(tmp %>% select(estimate), tmp %>% select(std.error))
      tibble(term = t, estimate = pool$q.mi[1], std.error = pool$se.mi[1]) 
    }) 
  }
  else{
    pool <- Amelia::mi.meld(imputed_model %>% select(estimate), imputed_model %>% select(std.error))
    tibble(estimate = pool$q.mi[1], std.error = pool$se.mi[1])
  }
}

na_frac <- function(x){
  mean(is.na(x))
}

simul_catr_error <- function(item_bank_file = "data/BAT_item_bank.csv", 
                             name = "BAT", 
                             size = 100, 
                             theta_mean = 0, 
                             theta_sd = 1,
                             num_items = 8,
                             seed = 666
                             ){
  if(num_items <= 0){
    check <- rnorm(size, theta_mean, theta_sd)
    if(any(is.na(check))){
      browser()
    }
    return(tibble(iter = 1:size, !!sym(sprintf("%s.score", name)) := rnorm(size, theta_mean, theta_sd), 
                  !!sym(sprintf("%s.error", name)) := 0))
    
  }
  set.seed(seed)
  alpha <- 0.05
  item_bank <- read.csv(item_bank_file, sep = ",")
  if(length(names(item_bank)) == 1){
    item_bank <- read.csv(item_bank_file, sep = ";")
  }
  #browser()
  if("difficulty_without_track_effect" %in% names(item_bank)){
    item_bank <- item_bank %>% rename(difficulty = difficulty_without_track_effect)
  }
  item_bank <- item_bank %>% select(discrimination, 
                                    difficulty, 
                                    guessing, 
                                    inattention) %>% as.matrix()
  start <- list(nrItems = 1, theta = 0, startSelect = "MFI")
  stop  <- list(rule = "length", thr = num_items, alpha = alpha)
  test  <- list(method = "BM", itemSelect = "bOpt", priorDist = "norm", priorPar = c(0, 1))
  final <- list(method = "WL", alpha = alpha)
  if(length(theta_mean) == 1){
    theta_mean <- rep(theta_mean, size)
  }
  if(length(theta_sd) == 1){
    theta_sd <- rep(theta_sd, size)
  }
  theta_mean[is.na(theta_mean)] <- mean(theta_mean, na.rm = T)
  theta_sd[is.na(theta_sd)] <- exp(mean(log(theta_sd), na.rm = T))
  stopifnot(length(theta_mean) == length(theta_sd))
  map_dfr(1:length(theta_mean), function(n){
    res <- catR::randomCAT(trueTheta = rnorm(1, theta_mean[n], theta_sd[n]), itemBank  = item_bank,
                           start = start, stop = stop, test = test, final = final)
    tibble(iter = n, score = res$thFinal, error = res$seFinal)
  }) %>% rename(!!sym(sprintf("%s.score", name)) := score, !!sym(sprintf("%s.error", name)) := error)
}

set_off_diag <- function(mat, off_diag_val = 0){
  
  if(!is.matrix(mat) || dim(mat)[1] != dim(mat)[2]){
    stop("*mat* must be quadratic matrix")
  }
  d <- dim(mat)[1L]
  res <- matrix(rep(off_diag_val, d^2), nrow = d)
  res[cbind(1L:d, 1L:d)] <- diag(mat)
  rownames(res) <- rownames(mat)
  colnames(res) <- colnames(mat)
  res
}

mult_diag <- function(mat, lambda = 0){
  if(!is.matrix(mat) || dim(mat)[1] != dim(mat)[2]){
    stop("*mat* must be quadratic matrix")
  }
  d <- dim(mat)[1]
  res <- mat
  res[cbind(1L:d, 1L:d)] <- lambda * diag(mat)
  res
}

simulate_with_corr <- function(data, size = nrow(data), sigma = NULL, var_scale = 1.0){

  stopifnot("MHE_class" %in% names(data))
  #browser()
  sum_stats <- data %>% 
    group_by(MHE_class) %>% 
    summarise_at(c("BAT.score", "MIQ.score", "MT.score", "MHE.score"), 
                 c(mean = mean, sd = sd), na.rm = T)
  
  map_dfr(unique(data$MHE_class), function(mhe_class){
    tmp <- sum_stats %>% filter(MHE_class == mhe_class)
    means <- c(tmp %>% pull(BAT.score_mean), 
               tmp %>% pull(MIQ.score_mean), 
               tmp %>% pull(MT.score_mean), 
               tmp %>% pull(MHE.score_mean))
    if(is.null(sigma)){
      sigma <- cov(data %>% filter(MHE_class == mhe_class) %>% 
                     select(BAT.score, MIQ.score, MT.score, MHE.score), use = "pairwise.complete.obs")
    }
    else  if(is.character(sigma)) {
      if(sigma == "decor"){
        sigma <- cov(data %>% filter(MHE_class == mhe_class) %>% 
                       select(BAT.score, MIQ.score, MT.score, MHE.score), use = "pairwise.complete.obs")
        sigma <- set_off_diag(sigma, off_diag_val = 0)
      }
      else if (sigma == "upcor"){
        sigma <- cov(data %>% filter(MHE_class == mhe_class) %>% 
                       select(BAT.score, MIQ.score, MT.score, MHE.score), use = "pairwise.complete.obs")
        sigma <- mult_diag(sigma, var_scale)
        if(rowSums(sigma)["MHE.score"] != 0){
          cov_test <- sigma %>% cov2cor()
        }
        else{
          cov_test <- sigma[1:3, 1:3] %>% cov2cor()
        }
        
        if(cov_test %>% is.na() %>% any() || (cov_test %>% (function(x) x > 1) %>% any())){
          stop(sprintf("Covariance matrix ill conditioned for scale factor %.3f", var_scale ))
        }
      }
    }
    else if(is.numeric(sigma)) {
        row.names(sigma) <- c("BAT.score", "MIQ.score", "MT.score", "MHE.score")
        colnames(sigma) <- c("BAT.score", "MIQ.score", "MT.score", "MHE.score")
      }
    else{
      stop("Sigma must be NULL, character or matrix")
    }

    n_effective <- round(size * nrow(data %>% filter(MHE_class == mhe_class))/nrow(data)) + 1
    MASS::mvrnorm(n = n_effective, mu = means, Sigma = sigma, empirical = T) %>% 
      as_tibble() %>% mutate(MHE_class = mhe_class)
  }) %>% `[`(1:size, )
}


cross_add_col <- function(df, col, col_name = NULL){
  if(is.null(col_name)){
    col_name <- deparse(substitute(col))
  }
  map_dfr(col, function(x){
    df %>% mutate(!!sym(col_name) := !!x)
  })
}

get_simu_def <- function(){
  method <- c("mvt", "mvt-decor", "mvt-upcor")
  MT.errors <- c(.76, .38, .318, 0)
  MIQ_items <- c(10, 20, 30, 0)
  BAT_items <- c(4, 8, 16, 0)
  with_na <- c(F, T)
  size <- c(250, 500, 1000, 2000, 4000)
  tmp <- tibble(MT.error = MT.errors, num_MIQ_items = MIQ_items, num_BAT_items = BAT_items, 
                error_level = c("large", "medium", "small", "none"), filename = "simu")
  cross_add_col(tmp, method) %>% cross_add_col(with_na) %>% cross_add_col(size) %>% mutate(id = 1:nrow(.))
}

simulate_data <- function(data = master_cross, 
                          size = nrow(data), 
                          seed = NULL, 
                          method = "mvt",
                          coefs = c(beta0 =  -1.1, beta_MT = 0.163, beta_MIQ = 0.355, beta_MHE = 0.016),
                          MT.error = .318,
                          num_MIQ_items = 8,
                          num_BAT_items = 20,
                          with_na = T, 
                          simu_params = NULL){
  #browser()
  if(!is.null(simu_params)){
    list2env(simu_params, environment())
  }
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(method == "faux"){
    tmp_data <- trim_data(data)
    tmp_data$MHE.score <- log(tmp_data$MHE.score - min(tmp_data$MHE.score) + 1)
    tmp <- faux::sim_df(tmp_data %>% select(all_of(measurement_vars)), 
                        n = size, 
                        id = "p_id", 
                        dv = "BAT.score")
    MIQ.true_score <- tmp$MIQ.score
    MT.true_score  <- tmp$MT.score
    MHE.true_score <- exp(tmp$MHE.score) 
    MIQ.theta_sd <- 0
  }
  else if(method == "mvt"){
    tmp <- simulate_with_corr(data, size)
    MT.true_score <- tmp$MT.score
    MHE.true_score <- tmp$MHE.score
    MIQ.true_score <- tmp$MIQ.score
    MIQ.theta_sd <- 0
  }
  else if(method == "mvt-decor"){
    tmp <- simulate_with_corr(data, size, sigma = "decor")
    MT.true_score <- tmp$MT.score
    MHE.true_score <- tmp$MHE.score
    MIQ.true_score <- tmp$MIQ.score
    MIQ.theta_sd <- 0
  }
  else if(method == "mvt-upcor"){
    tmp <- simulate_with_corr(data, size, sigma = "upcor", var_scale = .75)
    MT.true_score <- tmp$MT.score
    MHE.true_score <- tmp$MHE.score
    MIQ.true_score <- tmp$MIQ.score
    MIQ.theta_sd <- 0
  }
  else if (method == "manual"){
    MT.true_score  <- rnorm(size, mean = mean(data$MT.score, na.rm = T), sd = sd(data$MT.score, na.rm = T))
    MHE.true_score <- sample(data %>% trim_data(na.rm = T) %>% pull(MHE.score), size, replace = T)
    MIQ.true_score <- mean(data$MIQ.score, na.rm = T)  
    MIQ.theta_sd <- sd(data$MIQ.score, na.rm = T)  
  }
  else{
    stop(sprintf("Invalid method: %s", method))
  }
  #browser()
  tictoc::tic()
  messagef("Simulating MIQ values")
  miq_simul <- simul_catr_error(item_bank_file = "data/MIQ_item_bank.csv", 
                                name = "MIQ", 
                                size = size, 
                                theta_mean = MIQ.true_score, 
                                theta_sd = MIQ.theta_sd, 
                                num_items = num_MIQ_items) 

  tictoc::toc()

  MT.score <- (MT.true_score + rnorm(size, 0, MT.error))  %>% limiter(c(1, 7))
  
  BAT.true_score <- coefs["beta0"] + coefs["beta_MT"] * MT.true_score + coefs["beta_MIQ"] * miq_simul$MIQ.score + coefs["beta_MHE"] * MHE.true_score
  #fit_BAT <- lm(I(1/BAT.error) ~ (BAT.score), data = data) 
  tictoc::tic()
  messagef("Simulating BAT values")
  bat_simul <- simul_catr_error(item_bank_file = "data/BAT_item_bank.csv", 
                                size = size, 
                                theta_mean = BAT.true_score, 
                                theta_sd = 0, 
                                num_items = num_BAT_items)
  tictoc::toc()
  
  simul_data <- tibble(p_id = sprintf("S%d", 1:size), 
                  BAT.score = bat_simul$BAT.score,
                  BAT.true_score = BAT.true_score,
                  MT.score  = MT.score,
                  MT.error  = MT.error,
                  MIQ.score = miq_simul$MIQ.score,
                  MIQ.true_score = MIQ.true_score,
                  MIQ.error = miq_simul$MIQ.error,
                  MHE.score = MHE.true_score, 
                  BAT.error = bat_simul$BAT.error)
  #simul_data$BAT.error <- 1/predict(fit_BAT, simul_data) %>% 
  #  limiter(c(min(data$BAT.error, na.rm = T), max(data$BAT.error, na.rm = T)))
  #simul_data$BAT.score <- (BAT.true_score + rnorm(size, 0, simul_data$BAT.error)) %>% limiter(c(-4, 4))
  
  if(with_na){
    na_perc <-  map_dfr(all_vars, 
                      function(v) tibble(var = v, na_perc = na_frac(data[[var]]))) %>% 
      arrange(var)
    for(sv in score_vars){
      #browser()
      perc <- na_perc %>% filter(var == sv) %>% pull(na_perc)
      na_rows <- sample(1:nrow(simul_data), perc * nrow(simul_data))
      simul_data[na_rows, sv] <- NA
      simul_data[na_rows, gsub(".score", ".error", sv, fixed = T)] <- NA
    }
  }
  simul_data
}

setup_vars <- function(){
  assign("dep_var", "BAT.score", globalenv())
  assign("dep_var_error", "BAT.error", globalenv())
  assign("pred_primary", "MT.score", globalenv())
  assign("pred_primary_error", "MT.error", globalenv())
  assign("pred_secondary", "MIQ.score", globalenv())
  assign("pred_secondary_error", "MIQ.error", globalenv())
  assign("pred_tertiary", "MHE.score", globalenv())
  assign("dv_threshold", 1.596, globalenv())
  assign("measurement_vars", c(dep_var, pred_primary, pred_secondary, pred_tertiary), globalenv())
  assign("score_vars", c(dep_var, pred_primary, pred_secondary), globalenv())
  assign("error_vars", c(dep_var_error, pred_primary_error, pred_secondary_error), globalenv())
  assign("simex_vars", c(pred_primary, pred_secondary), globalenv())
  assign("simex_error", c(pred_primary_error, pred_secondary_error), globalenv())
  assign("all_vars", c(measurement_vars, error_vars), globalenv())
}

setup_models <- function(){
  predictors <- paste(list(pred_primary, pred_secondary, pred_tertiary), collapse = " + ")
  assign("lm_model", as.formula(sprintf("%s ~ %s", dep_var, predictors)), globalenv())
  assign("brms_model_se", as.formula(sprintf("%s | resp_se(%s) ~ %s", dep_var, dep_var_error, predictors)), globalenv())
  assign("brms_model_ses", as.formula(sprintf("%s | resp_se(%s, sigma = TRUE) ~ %s", dep_var, dep_var_error, predictors)), globalenv())
}

setup_workspace_me <- function(reload_data = T){
  messagef("Setting up vars...")
  setup_vars()
  messagef("Setting up models...")  
  setup_models()
  
  if(reload_data){
    aws.set_credentials_if_unset()
    messagef("Reloading and preparing data...")
    setup_workspace()
    master_cross <- get_cross_section_version(master)
    master_cross <- master_cross %>% 
      rename(MHE.score = MHE.general_score, 
             MT.score = GMS.musical_training) %>% 
      mutate(MT.error = .318) %>% 
      trim_data(na.rm = F) %>% 
      filter(Reduce(`+`, lapply(., is.na)) != ncol(.) - 2) 

    master_cross[is.na(master_cross$MT.score), ]$MT.error <- NA
    data <- mice::mice(master_cross %>% select(BAT.score, MIQ.score, MT.score, MHE.score), m = 1, method ="pmm", printFlag = F) %>% complete(1)
    #browser()
    #data[is.na(data$MHE.score),]$MHE.score <- imp$imp$MHE.score[[1]]
    MHE <- data$MHE.score
    MHE_GMM <- mclust::Mclust(log(MHE - min(MHE) + .1) + rnorm(length(MHE), 0, .01), G = 3, na.rm = T)
    master_cross$MHE_class <- MHE_GMM$classification 
  }
  else{
    messagef("Loading data...")
    master_cross <- readRDS("data/master_cross.rds")
  }
  messagef("Done.")
  
  assign("master_cross", master_cross, globalenv())
}

trim_data <- function(data, na.rm = T){
  if(na.rm == T){
    data <- data %>% filter(!is.na(!!sym(dep_var)),
                    !is.na(!!sym(pred_primary)),
                    !is.na(!!sym(pred_primary_error)),
                    !is.na(!!sym(pred_secondary)),
                    !is.na(!!sym(pred_secondary_error)),
                    !is.na(!!sym(pred_tertiary))) 
    
  }
  data %>% select(p_id, 
                  all_of(dep_var), 
                  all_of(pred_primary), 
                  all_of(pred_primary_error), 
                  all_of(dep_var_error), 
                  all_of(pred_secondary), 
                  all_of(pred_secondary_error), 
                  all_of(pred_tertiary))
}
add_confint_tidy <- function(lm_fit, level = .95){
  tidy <- broom::tidy(lm_fit)
  cis <- map_dfr(attr(lm_fit$coefficients, "names"), function(x){
    confint(lm_fit, x, level) %>% 
      as_tibble() %>% 
      rename(low_ci = 1, up_ci = 2) %>% 
      mutate(term = x)
    })
  tidy %>% left_join(cis, by = "term")
}

lm_wrapper <- function(data, broom_FUN =  broom::tidy, use_simex = F, ...){
  
  data <- trim_data(data, na.rm = T)
  fun_name <- deparse(substitute(broom_FUN))
  
  if(use_simex){
    #browser()
    fit_raw <- lm(lm_model, x = T, data = data)  
    arguments <- list(...)
    arg_names <- names(arguments)
    if(length(intersect(c("simex_var", "simex_error"), names(arguments))) != 2){
      stop("If use_simex = T you must specify 'simex_var' and 'simex_error'")
    }   
    simex_var <- arguments[["simex_var"]]
    simex_error <- arguments[["simex_error"]]
    if(sum(data[simex_error] == 0)){
      fit <- NULL
    }
    else{
      fit<- simex::simex(fit_raw, measurement.error = data[simex_error], 
                            SIMEXvariable = simex_var, asymptotic = F) 
    }
    if(!is.null(fit) && !is.null(broom_FUN)){
      fit  <- summary(fit)
      fit <- fit$coefficients$jackknife %>% as.data.frame()
      fit <- rownames_to_column(fit)
      names(fit) <- c("term", "estimate", "std.error", "statistic", "p.value")
      fit <- fit %>% mutate(low_ci = estimate - 1.96*std.error, up_ci = estimate + 1.96*std.error)
    }
  }
  else{
    if(fun_name == "broom::glance"){
      fit <- lm(lm_model, data = data) %>% broom_FUN()
    }
    else{
      fit <- lm(lm_model, data = data) %>% add_confint_tidy()
    }
  }
  #browser()
  fit
}

test_overimputation <- function(data, m = 5){
  set.seed(666)
  #data <- trim_data(data)
  #data <- data %>% sample_n(30)
  messagef("***Calculating overimputation")
  #browser()
  data <- data %>% select(-ends_with("true_score")) 

  overimp_mat <- tibble(row = integer(0), column = integer(0))
  priors <- tibble(row = integer(0), column = integer(0), x = double(0), y = double(0))
  for(iv in score_vars){
    var_column <- which(names(data) == iv)
    good_rows <- which(!is.na(data[[iv]])) 
    tmp_mat <- tibble(row = good_rows, column = var_column)
    overimp_mat <- bind_rows(overimp_mat, tmp_mat)
    iv_error <- gsub(".score", ".error", iv, fixed = T)
    priors <- bind_rows(priors, 
                        bind_cols(tmp_mat, 
                                  data[good_rows,] %>% 
                                    select(x = all_of(iv), y = all_of(iv_error))))
    
  }
  #browser()
  #data <- data %>% select(-dep_var)
  id_vars <- c("p_id", "MT.error")
  if(sd(data$BAT.error, na.rm = T) <= 0.01){
    id_vars <- c(id_vars, "BAT.error")
  }
  if(sd(data$MIQ.error, na.rm = T) <=   0.01){
    id_vars <- c(id_vars, "MIQ.error")
  }
  if(sd(priors$y, na.rm = T) <=   0.01){
    priors$y <- rnorm(nrow(priors), 0, 0.001)
  }
  
  imputed <- Amelia::amelia(data %>% as.data.frame(),
                    m = m,
                    idvars = id_vars, 
                    overimp = as.matrix(overimp_mat),
                    priors = as.matrix(priors), 
                    p2s = 0
                    )
  #browser()
  assign("imputed", imputed, globalenv())
  imputed_lm(imputed, lm_model) %>% 
    mutate(type = "overimputation") %>% 
    select(term, estimate, std.error, statistic, p.value, low_ci, up_ci, type)
}




test_pmm_imputation <- function(data, m = 5){

  set.seed(666)
  
  messagef("***Calculating pmm imputation")
  
  #browser()
  data <- data %>% select(-ends_with("true_score"))
  imp <- mice::mice(data %>% select(-all_of(error_vars)), m = m, printFlag = F, method = "pmm")
  
  fit <- 
    map(1:m, function(i) {
      lm(lm_model, data = mice::complete(imp, i)) 
    })
  
  summary(fit  %>% mice::pool()) %>% 
    mutate(low_ci = estimate - 1.96 * std.error, up_ci = estimate + 1.96 * std.error, type = "pmm_imputation")
}

test_simex <- function(data){
  set.seed(666)
  fit <- lm_wrapper(data, use_simex = T, simex_var = simex_vars, simex_error = simex_error) 
  if(!is.null(fit)) fit <- fit %>%  mutate(type = "simex_raw")
  fit
}
get_model_summary <- function(data, model, type, cutoff = NA, use_weights = F, ...){
  data <- trim_data(data)
  if(use_weights){
    #data$weights <- 1/data[[dep_var_error]]^2
    if(sum(data[, error_vars]) == 0){
      return(NULL)
    }
    data$weights <- 1/sqrt(rowSums(data[, error_vars]^2)) 
    fit <- lm(model, data = data, weights = I(weights)) 
  }
  else{
    fit <- lm(model, data = data) 
  }
  fit %>%  add_confint_tidy() %>% mutate(type = type)  
  
}

test_all <- function(data, m = 5){
  fit_overimp <- test_overimputation(data, m = m)
  fit_pmm_imp <- test_pmm_imputation(data, m = m)

  fit_simex <- test_simex(data)
  fit_raw <- get_model_summary(data, lm_model, type = "raw")
  fit_raw_w <- get_model_summary(data, lm_model, type = "weighted", use_weights = T)
  fit_true <- tibble(term = c("(Intercept)", "MT.score", "MIQ.score", "MHE.score"), 
                     estimate = c(-1.1, .166, .355, 0.016),
                     low_ci = c(-1.1, .166, .355, 0.016),
                     up_ci = c(-1.1, .166, .355, 0.016),
                     type = "true")
  pool <- 
    bind_rows(
      fit_true,
      fit_raw,
      fit_raw_w,
      fit_overimp, 
      fit_pmm_imp, 
      fit_simex) 
  pool
}

in_interval <- function(x, low, up){
  any(is.na(c(low, up))) || ((x>= low) && (x <= up))
}

get_relative_stats <- function(pool){
  pool %>% 
    filter(type == "true") %>% 
    rename(true_beta = estimate) %>% 
    distinct(term, true_beta) %>% 
    left_join(pool %>% select(term, estimate, type, low_ci, up_ci), by = "term") %>% 
    mutate(attenuation = estimate/true_beta, 
           abs_diff = abs(estimate-true_beta), 
           abs_rel_diff = abs(estimate-true_beta)/true_beta,
           rel_diff = (estimate-true_beta)/true_beta,
           in_ci = in_interval(true_beta, low_ci, up_ci),
           ci_coverage = as.numeric(in_ci)) %>%  
    arrange(term, type) %>% select(term, type, everything())
}

test_simulations <- function(data, n_simul = 30, imp_m = 5, simu_params = NULL, label = "DEFAULT"){
  if(is.null(simu_params)){
    simu_params <- get_simu_def()[1,]
  }
  #browser()

  raw <-  
    map_dfr(1:n_simul, function(n){
      deco_messagef("%s: Simulating data set #%d/%d",  label, n, n_simul)
      simu <- simulate_data(data, simu_params = simu_params)
      simu %>% mutate(iter = n)
  })
  
  pool <- 
    map_dfr(1:n_simul, function(n){
      deco_messagef("%s: Testing data set #%d/%d", label, n,  n_simul)
      test_all(raw %>% filter(iter == n) %>% select(-iter), m = imp_m) %>% mutate(iter = n)
  })
  #browser()
  stats <- 
    map_dfr(1:n_simul, function(n){
      #browser()
      deco_messagef("%s: Collecting stats for set #%d/%d", label, n, n_simul)
      rel_stats <- get_relative_stats(pool %>% filter(iter == n, term != "(Intercept)")) %>% 
        mutate(iter = n) %>% select(-low_ci, -up_ci)
      sum_stats <- 
        rel_stats %>% 
          group_by(type) %>% 
          summarise(abs_diff = mean(abs_diff), 
                    rel_diff = mean(rel_diff),
                    abs_rel_diff = mean(abs_rel_diff),
                    attenuation = mean(attenuation),
                    ci_coverage = mean(in_ci, na.rm = T), .groups = "drop") %>% 
          ungroup() %>% 
          mutate(iter = n, term = "sum_stats")
      bind_rows(rel_stats, sum_stats)
      
    })
  list(raw = raw, 
       pool = pool, 
       stats = stats, 
       params = simu_params %>% mutate(n_simul = n_simul, imp_m = imp_m), 
       version = version)
}

test_all_simulations <- function(data, n_simul, imp_m = 5, 
                                 simu_def = NULL, 
                                 label = "simu", 
                                 out_dir = "data/simulations"){
  if(is.null(simu_def)){
    simu_def <- get_simu_def()[1,]
  }
  save_data <-  is.character(out_dir) && nchar(out_dir) > 0
  ret <- list()
  for(r in 1:nrow(simu_def)){
    tictoc::tic()
    params <- sprintf("%s: %s", simu_def[r,] %>% names(), simu_def[r,] %>% as.list()) %>% paste(collapse ="\n")
    deco_messagef("Running %s", sprintf("%s/%d", label, r))
    tmp <- test_simulations(data = data, n_simul = n_simul, imp_m = imp_m, 
                            simu_params = simu_def[r, ], label = sprintf("%s/param_id=%d", label, simu_def[r,]$id))
    tictoc::toc()
    #browser()
    if(save_data){
      filename <- file.path(out_dir,
                             sprintf("%s_meth=%s_sz=%d_el=%s_n=%d_m=%d_NA=%d.rds", 
                                     simu_def[r,]$filename, 
                                     simu_def[r, ]$method, 
                                     simu_def[r, ]$size, 
                                     simu_def[r, ]$error_level,
                                     n_simul, 
                                     imp_m,
                                     as.integer(simu_def[r, ]$with_na)))
      filename <- file.path(out_dir, sprintf("%s_id=%02d_n=%d.rds", label, simu_def[r,]$id, n_simul))
      saveRDS(tmp, filename)
      messagef("Uploading results to AWS...%s", filename)
      #browser()
      aws.s3::s3saveRDS(tmp, 
                        object = sprintf("%s/%s", aws3_dir, filename),
                        bucket = "longgold.gold-msi.org",
                        opts = list(region = "eu-west-1"))
      
    }
    ret[[r]] <- tmp
  }

  #raw <- map_dfr(ret, function(x) x$raw %>% bind_cols(x$params) %>% bind_cols(x$version))
  pool <- map_dfr(ret, function(x) x$pool %>% bind_cols(x$params) %>% mutate(version = version))
  stats <- map_dfr(ret, function(x) x$stats %>% bind_cols(x$params) %>% mutate(version = version))
  simu_data = list(pool = pool, stats = stats)
  deco_messagef("%s: Done!", label)
  
  if(save_data){
    filename <- file.path(out_dir, sprintf("%s_v=%s_%s.rds", label, version, digest::sha1(lubridate::now())))
    deco_messagef("%s: Saving all summaries to: '%s'", label, filename)
    saveRDS(simu_data, file = filename) 
    #browser()
    aws.s3::s3saveRDS(simu_data, 
                      object = sprintf("%s/%s", aws3_dir, filename),
                      bucket = "longgold.gold-msi.org",
                      opts = list(region = "eu-west-1"))
  }
  else{
    simu_data
  }
}
