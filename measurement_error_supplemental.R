
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
  f <- as.formula("BAT.score ~ MT.score + MIQ.score")
  raw_lm <- lm(f, data = data) %>% broom::tidy()
  #browser()
  brat_wrapper <- function(data){
    lm(f, data = data) %>% 
      broom::tidy() 
  }
  lm_brat_pred <- bootstrapper(data, 
                               FUN = brat_wrapper, 
                               vars = c("MT"), 
                               size = size)
  
  brat_pred <- pool_boots_lm(lm_brat_pred) %>% mutate(model = "kf_impute_pred")
  
  lm_brat_dep <- bootstrapper(data, 
                              FUN = brat_wrapper, 
                              vars = c("BAT"), 
                              size = size)
  
  brat_dep <- pool_boots_lm(lm_brat_dep) %>% mutate(model = "kf_impute_dep")
  
  lm_brat <- bootstrapper(data, 
                          FUN = brat_wrapper, 
                          vars = c("BAT", "MT"), 
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