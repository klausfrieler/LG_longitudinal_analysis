library(tidyverse)

version <- "0.1.0"

#This hack avoids including library(mclust), due to bug in mclust
mclustBIC <- mclust::mclustBIC

source("setup.R")
#source("scripts/IRT_bootstrapper.R")

required_packages <- c("aws.s3", "lubridate", 
                       "brms", "mice", "miceadds", "mclust", 
                       "catR", "Amelia", "broom", "MASS", 
                       "faux", "broom.mixed", "simex")

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
    return(tibble(iter = 1:size, !!sym(sprintf("%s.score", name)) := theta_mean, 
                  !!sym(sprintf("%s.error", name)) := rnorm(size, 0, .001)))
    
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
                 c(mean = mean,sd = sd), na.rm = T)
  
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

    n_effective <- round(size * nrow(data %>% filter(MHE_class == mhe_class))/nrow(data))
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
  MT.errors <- c(.76, .38, .318)
  MIQ_items <- c(10, 20, 30)
  BAT_items <- c(4, 8, 16)
  with_na <- c(F, T)
  size <- c(250, 500, 1000, 2000, 4000, 8000)
  tmp <- tibble(MT.error = MT.errors, num_MIQ_items = MIQ_items, num_BAT_items = BAT_items, 
                error_level = c("large", "medium", "small"), filename = "simu")
  cross_add_col(tmp, method) %>% cross_add_col(with_na) %>% cross_add_col(size)
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
    data <- mice::mice(master_cross %>% select(BAT.score, MIQ.score, MT.score, MHE.score), m = 1, method ="pmm") %>% complete(1)
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
  tidy %>% left_join(cis)
}

lm_wrapper <- function(data, broom_FUN =  broom::tidy, use_simex = F, ...){
  data <- trim_data(data, na.rm = T)

  fun_name <- deparse(substitute(broom_FUN))
  if(use_simex){
    fit_raw <- lm(lm_model, x = T, data = data)  
    arguments <- list(...)
    arg_names <- names(arguments)
    if(length(intersect(c("simex_var", "simex_error"), names(arguments))) != 2){
      stop("If use_simex = T you must specify 'simex_var' and 'simex_error'")
    }   
    simex_var <- arguments[["simex_var"]]
    simex_error <- arguments[["simex_error"]]
    fit<- simex::simex(fit_raw, measurement.error = data[simex_error], 
                            SIMEXvariable = simex_var, asymptotic = F) 
    if(!is.null(broom_FUN)){
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
  fit
}

get_model_summary <- function(data, model, type, cutoff = NA, use_weights = F, ...){
  data <- trim_data(data)
  if(use_weights){
    #data$weights <- 1/data[[dep_var_error]]^2
    data$weights <- 1/sqrt(rowSums(data[, error_vars]^2)) 
    fit <- lm(model, data = data, weights = I(weights)) 
  }
  else{
    fit <- lm(model, data = data) 
  }
  fit %>%  broom::tidy() %>% mutate(type = type)
 
}

test_cutoff <- function(data, model, error_var = dep_var_error, dv_threshold = NULL){
  set.seed(666)
  #data <- trim_data(data, na.rm = F)
  messagef("***Calculating cutoffs")
  
  l <- data %>% nrow()
  fit_raw <- get_model_summary(data, model, type = "raw")
  fit_raw_w <- get_model_summary(data, model, type = "weighted", use_weights = T)
  if(is.null(dv_threshold)){
    min_t <- round(min(data[[dep_var_error]], na.rm = T), 1)
    max_t <- round(max(data[[dep_var_error]], na.rm = T), 1)
    dv_threshold <- seq(min_t + .01, max_t, .1)
  }
  map_dfr(dv_threshold, function(cutoff){
    messagef("Testing cutoff %.3f", cutoff)
    #browser()
    data_cutoff <- data %>%
      filter(!!sym(error_var) <= cutoff) 
    l_cutoff <- nrow(data_cutoff)
    data_sampling <- data %>% sample_frac(l_cutoff/l, replace = F)
    if(nrow(data %>% filter(!is.na(!!sym(dep_var)))) == 0){
      return(NULL)
    }
    fit_sampling <- get_model_summary(data = data_sampling, model, type = "sampling", cutoff = cutoff) 
    fit_cutoff <- get_model_summary(data = data_cutoff, model, type = "cutoff", cutoff = cutoff) 
    bind_rows(fit_sampling, fit_cutoff)
  }) %>% bind_rows(fit_raw, fit_raw_w)
}

test_brms <- function(data, iter = 2500){
  set.seed(666)
  data <- trim_data(data) 
  messagef("***Calculating brms models")
  
  fit <- brms::brm(lm_model, 
                   data = data, 
                   cores = 4, 
                   iter = iter)
  fit_se <- brms::brm(brms_model_se, data = data, cores = 4, iter = iter)
  fit_ses <- brms::brm(brms_model_ses, data = data, cores = 4, iter = iter)
  #browser()
  
  bind_rows(fit %>% broom.mixed::tidyMCMC() %>% mutate(type = "brms_raw"), 
            fit_se %>% broom.mixed::tidyMCMC() %>% mutate(type = "brms_se"),
            fit_ses %>% broom.mixed::tidyMCMC() %>% mutate(type = "brms_ses")) %>% 
    mutate(term = gsub("b_", "", term)) %>% 
    mutate(term = gsub("Intercept", "(Intercept)", term)) %>% 
    filter(term != "sigma")
    
}

se <- function(x,...){
  sd(x, ...)/sqrt(length(x))
}


test_imputation_KF <- function(data, size = 30){
  set.seed(666)
  messagef("***Calculating imputation KF")
  data <- data %>% trim_data(na.rm = F)
  
  boot_lm <- bootstrapper(data, FUN = lm_wrapper, 
                          score_vars = score_vars, 
                          error_vars = error_vars, 
                          size = size, 
                          broom_FUN = broom::tidy)

  boot_lm <- quick_meld(boot_lm, by_term = T)

  raw_lm_simex <- data %>% 
    lm_wrapper(broom_FUN = broom::tidy, use_simex = T, simex_var = simex_vars,  simex_error = simex_error) 
  
  boot_lm_simex <- bootstrapper(data, FUN = lm_wrapper, 
                          score_vars = dep_var, 
                          error_vars = dep_var_error, 
                          use_simex = T, 
                          simex_var = simex_vars, 
                          simex_error = simex_error,
                          size = size, 
                          broom_FUN = broom::tidy) 

  boot_lm_simex <- quick_meld(boot_lm_simex, by_term = T)
  
  bind_rows(boot_lm %>% mutate(type = "boot_lm"), 
            #raw_lm %>% mutate(type = "raw_lm"),
            boot_lm_simex %>% mutate(type = "boot_lm_simex"), 
            raw_lm_simex %>% mutate(type = "raw_lm_simex")
            ) %>% select(term, everything())
  
}

test_overimputation <- function(data, m = 5){
  set.seed(666)
  #data <- trim_data(data)
  #data <- data %>% sample_n(30)
  messagef("***Calculating overimputation")
  #browser()
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
  if(sd(data$BAT.error) <= 0.01){
    id_vars <- c(id_vars, "BAT.error")
    
  }
  if(sd(data$MIQ.error) <=   0.01){
    id_vars <- c(id_vars, "MIQ.error")
    
  }
  imputed <- Amelia::amelia(data %>% as.data.frame(),
                    m = m,
                    idvars = id_vars, 
                    overimp = as.matrix(overimp_mat),
                    priors = as.matrix(priors), 
                    p2s = 0,
                    ncpus = 8,
                    )
  #browser()
  assign("imputed", imputed, globalenv())
  imputed_lm(imputed, lm_model) %>% 
    mutate(type = "overimputation") %>% 
    select(term, estimate, std.error, statistic, p.value, low_ci, up_ci, type)
}


test_pv_imputation <- function(data, m = 30){

  set.seed(666)
  
  #data <- trim_data(data, na.rm = T)
  #data <- data %>% select(p_id, all_of(dep_var), all_of(pred_secondary))
  messagef("***Calculating plausible value imputation")
  
  data[["p_id"]] <- NULL

  params <- list("dep_var" = list("M" = data[[dep_var]], "SE" = data[[dep_var_error]]),
                 "pred_primary" = list("M" = data[[pred_primary]], "SE" = data[[pred_primary_error]]),
                 "pred_secondary" = list("M" = data[[pred_secondary]], "SE" = data[[pred_secondary_error]]))
  names(params) <- c(dep_var, pred_primary, pred_secondary)
  data[[dep_var_error]] <- NULL
  data[[pred_primary_error]] <- NULL
  data[[pred_secondary_error]] <- NULL
  cn <- colnames(data)
  V <- length(cn)
  method <- rep("", V)
  names(method) <- cn
  method[dep_var] <- "plausible.values"
  method[pred_primary] <- "plausible.values"
  method[pred_secondary] <- "plausible.values"
  
  browser()
  imp <- mice::mice(data, 
                    maxit = 1, 
                    m = m, 
                    allow.na = TRUE, 
                    method = method, 
                    printFlag = FALSE,
                    scale.values = params)

  
  #fit <- mice::pool(with(imp, lm(eval(lm_model))))
  #lm_model_loc <- as.formula(sprintf("%s~%s", dep_var, pred_secondary))
  fit <- 
    purrr:::map(1:m, function(i) {
      lm(lm_model, data = mice::complete(imp, i))
    }) %>% mice::pool()
  #fit <- fit$pooled %>% filter(term == pred_primary)
  summary(fit) %>% mutate(type = "pv_imputation")
}

test_pmm_imputation <- function(data, m = 30){

  set.seed(666)
  
  messagef("***Calculating plausible value imputation")
  
  #browser()
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
  lm_wrapper(data, use_simex = T, simex_var = simex_vars, simex_error = simex_error) %>% 
    mutate(type = "simex_raw")

}

test_all <- function(data, m = 5){
  #fit_cutoff <- test_cutoff(data, lm_model, dv_threshold = dv_threshold)
  #fit_brms <- test_brms(master_cross)
  #fit_brms <- NULL
  #fit_imp_kf <- test_imputation_KF(data) 
  fit_overimp <- test_overimputation(data, m = m)
  fit_pmm_imp <- test_pmm_imputation(data, m = m)
  fit_simex <- test_simex(data)
  fit_true <- tibble(term = c("(Intercept)", "MT.score", "MIQ.score", "MHE.score"), 
                     estimate = c(-1.1, .166, .355, 0.016),
                     low_ci = c(-1.1, .166, .355, 0.016),
                     up_ci = c(-1.1, .166, .355, 0.016),
                     type = "true")
  pool <- 
    bind_rows(
      #fit_brms, 
      #fit_cutoff, 
      #fit_imp_kf, 
      fit_true,
      fit_overimp, 
      fit_pmm_imp, 
      fit_simex) 
  pool
}

get_relative_stats <- function(pool){
  pool %>% 
    filter(type == "true") %>% 
    rename(true_beta = estimate) %>% 
    select(term, true_beta) %>% 
    left_join(pool %>% select(term, estimate, type), by = "term") %>% 
    mutate(attenuation = estimate/true_beta, 
           abs_diff = abs(estimate-true_beta), 
           abs_rel_diff = abs(estimate-true_beta)/true_beta,
           rel_diff = (estimate-true_beta)/true_beta) %>%  
    arrange(term, type) %>% select(term, type, everything())
  
}

test_simulations <- function(data, n_simul = 30, imp_m = 5, simu_params = NULL){
  if(is.null(simu_params)){
    simu_params <- get_simu_def()[1,]
  }
  #browser()
  raw <-  
    map_dfr(1:n_simul, function(n){
      simu <- simulate_data(data, simu_params = simu_params)
      simu %>% mutate(iter = n)
  })
  
  pool <- 
    map_dfr(1:n_simul, function(n){
      #browser()
      test_all(raw %>% filter(iter == n) %>% select(-iter), m = imp_m) %>% mutate(iter = n)
  })
  #browser()
  stats <- 
    map_df(1:n_simul, function(n){
      #browser()
      rel_stats <- get_relative_stats(pool %>% filter(iter == n, term != "(Intercept)")) %>% mutate(iter = n)
      sum_stats <- 
        rel_stats %>% 
          group_by(type) %>% 
          summarise(abs_diff = mean(abs_diff), 
                    rel_diff = mean(rel_diff),
                    abs_rel_diff = mean(abs_rel_diff),
                    attenuation = mean(attenuation)) %>% 
          ungroup() %>% 
          mutate(iter = n, term = "sum_stats")
      bind_rows(rel_stats, sum_stats)
      
    })
  list(raw = raw, pool = pool, stats = stats, params = simu_params %>% mutate(n_simul = n_simul, imp_m = imp_m), version = version)
}

test_all_simulations <- function(data, n_simul, imp_m = 5, simu_def = NULL, out_dir = "data/simulations"){
  if(is.null(simu_def)){
    simu_def <- get_simu_def()[1,]
  }
  ret <- list()
  for(r in 1:nrow(simu_def)){
    tictoc::tic()
    tmp <- test_simulations(data = data, n_simul = n_simul, imp_m = imp_m, simu_params = simu_def[r, ])
    tictoc::toc()
    file_name <- file.path(out_dir,
                           sprintf("%s_meth=%s_sz=%d_el=%s_n=%d_m=%d.rds", 
                                   simu_def[r,]$filename, 
                                   simu_def[r, ]$method, 
                                   simu_def[r, ]$size, 
                                   simu_def[r, ]$error_level,
                                   n_simul, 
                                   imp_m))
    saveRDS(tmp, file_name)
    ret[[r]] <- tmp
  }
  #raw <- map_dfr(ret, function(x) x$raw %>% bind_cols(x$params) %>% bind_cols(x$version))
  pool <- map_dfr(ret, function(x) x$pool %>% bind_cols(x$params) %>% mutate(version = version))
  stats <- map_dfr(ret, function(x) x$stats %>% bind_cols(x$params) %>% mutate(version = version))
  simu_data = list(pool = pool, stats = stats)
  save(simu_data, file = file.path(out_dir, sprintf("simu_%s.rda", version)))  
}

vary_tertiaries <- function(){
  orig_tertiary <- pred_tertiary
  tertiaries <- c("CCM", "MHE.general_score", "SES.educational_degree", "GMS.active_engagement")
  ret <- 
    map_dfr(tertiaries, function(x){
      assign("pred_tertiary", x, globalenv())
      test_all() %>% mutate(model = x)
      }) %>% 
    mutate(type = fct_reorder(factor(type), rel_diff, mean))
  assign("pred_tertiary", orig_tertiary, globalenv())
  ret
}

pool_boots_lm <- function(lm_model, dep_var  = "RAT.score"){
  #browser()
  pred_pool <- Amelia::mi.meld(lm_model %>% filter(term == dep_var) %>% select(estimate), 
                               lm_model %>% filter(term == dep_var) %>% select(std.error))
  intercept_pool <- Amelia::mi.meld(lm_model %>% filter(term == "(Intercept)") %>% select(estimate), 
                                    lm_model %>% filter(term == "(Intercept)") %>% select(std.error)) 
  #browser()
  bind_rows(
    tibble(term = "(Intercept)", 
           estimate = intercept_pool$q.mi[1], 
           std.error = intercept_pool$se.mi[1]),
    tibble(term = "(RAT.score)", 
           estimate = pred_pool$q.mi[1], 
           std.error = pred_pool$se.mi[1]))
  
}

test_IRT_kf_imputation <- function(data = master_cross, size = 100){
  f <- as.formula("BAT.score ~ RAT.score + MIQ.score")
  raw_lm <- lm(f, data = data) %>% broom::tidy()
  #browser()
  brat_wrapper <- function(data){
    lm(f, data = data) %>% 
      broom::tidy() 
  }
  lm_brat_pred <- bootstrapper(data, 
                       FUN = brat_wrapper, 
                       vars = c("RAT"), 
                       size = size)
  
  brat_pred <- pool_boots_lm(lm_brat_pred) %>% mutate(model = "kf_impute_pred")

  lm_brat_dep <- bootstrapper(data, 
                               FUN = brat_wrapper, 
                               vars = c("BAT"), 
                               size = size)
  
  brat_dep <- pool_boots_lm(lm_brat_dep) %>% mutate(model = "kf_impute_dep")
  
  lm_brat <- bootstrapper(data, 
                              FUN = brat_wrapper, 
                              vars = c("BAT", "RAT"), 
                              size = size)
  
  brat_pred_dep <- pool_boots_lm(lm_brat) %>% mutate(model = "kf_impute_pred_dep")
  #browser()
  bind_rows(raw_lm %>% select(term, estimate, std.error) %>% mutate(model = "raw"),
            brat_pred,
            brat_dep,
            brat_pred_dep)
  
}
test_xyz <- function(batch_size = 100, 
                     size = 1000, 
                     m = 0, 
                     s = 1){
  map_dfr(1:size, function(j){
    x <- rnorm(batch_size, 0, s)
    m <- mean(x)
    se <- se(x)
    y <- rnorm(batch_size, m, se)
    tibble(iter = j, 
           true_mean = m , 
           true_sd = s, 
           batch = batch_size, 
           size = size, 
           sim_mean = mean(y), 
           sim_sd = sd(y),
           d_mean = true_mean - sim_mean,
           d_sd = true_sd - sim_sd) 
  })
  
}