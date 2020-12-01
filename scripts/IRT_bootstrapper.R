source("setup.R")
#setup_workspace()

sem_wrapper <- function(data, ...){
  model_music <-
    '# measurement model
  model_music =~ MDI.score + MPT.score + BAT.score'
  lavaan::sem(model_music, data) %>% lavaan::standardizedsolution() 
}

lm_wrapper <- function(data, broom_FUN =  broom::glance, ...){
  lm(academic.overall ~ MIQ.score + RAT.score, data = data) %>% broom_FUN()   
}

cor_test_wrapper <- function(data){
  cor.test( ~ MIQ.score + RAT.score, data = data) %>% broom::tidy()
}
stat_desc_wrapper <- function(var, FUN.stat = "mean", ...){
  function(data){
    #browser()
    stat_name <- sprintf("%s_%s", FUN.stat, var)
    tibble(!!sym(stat_name) := eval(parse(text = FUN.stat))(data %>% pull(!!sym(var)), ...))      
    
  }
}

id_wrapper <- function(var){
  function(data){
    data %>% select(!!sym(var))
  }
}  
bootstrapper <- function(data, FUN, vars = NULL, score_vars = NULL, error_vars = NULL, size = 100, limits = c(-4,4), ...){
  if(!is.null(vars)){
    score_vars <- sprintf("%s.score", vars)
    error_vars <- sprintf("%s.error", vars)
  }
  else{
    stopifnot(!is.null(score_vars) & !is.null(error_vars))
  }
  var_pool <- 
    map(1:length(score_vars), function(i){
    list(scores = data %>% pull(score_vars[i]), errors = data %>% pull(error_vars[i]))
  })
  #browser()
  map_dfr(1:size, function(i){
    #browser()
    #messagef("Bootstrap iteration %d", i)
    map(1:length(var_pool), function(j){
      #messagef("Doing %s", score_vars[[j]])
      x_boot <-
        map2_dbl(var_pool[[j]]$scores, var_pool[[j]]$errors, function(m, s){
          suppressWarnings(rnorm(1, m, s)) %>% limiter(limits = limits)
        })
      data[[score_vars[j]]] <<- x_boot
    })
    FUN(data = data, ...) %>% mutate(type = "boot", iteration = i)
    #cor.test(x_boot, y_boot) %>% broom::tidy() %>% select(estimate, p.value) %>% mutate(type = "boot")
  })    
}

distribution_plot <- function(data, var, raw_value, group_var = "", title = ""){
  #browser()
  q <- data %>% ggplot(aes(x = !!sym(var), y =..count..)) + geom_histogram()
  if(nchar(group_var) != 0){
    line_data <- bind_cols(
      tibble(
      raw_value = raw_value),
      data %>% 
        group_by(!!sym(group_var)) %>% 
        summarise(mean_value = mean(!!sym(var), na.rm =  T)) %>% 
        ungroup() 
      ) 
   q <- q + geom_vline(data = line_data, aes(group = rhs, xintercept = raw_value), colour = "green")
   q <- q + geom_vline(data = line_data, aes(group = rhs, xintercept = mean_value), colour = "red")
   q <- q  + facet_wrap(as.formula(sprintf("~%s", group_var)))
   
  }
  else{
    mean_value <- mean(data[[var]], na.rm = T)
    q <- q + geom_vline(xintercept = raw_value, colour = "green")
    q <- q + geom_vline(xintercept = mean_value, colour = "red")
  }

  q <- q  + theme_bw()
  if(is.null(title) || length(title) == 0 || nchar(title) == 0){
    title <- var
  }
  q <- q + labs(x = var, title = title)
 
  q
}

demo_bootstrapper <- function(data = master %>% filter(test_year == 2019), size = 100){

  #1. check correlation
  messagef("Running correlation bootstrap with %d iterations.", size)
  raw_cor <- data %>% cor_test_wrapper() %>% pull(estimate)
  boot_cor <- bootstrapper(data, FUN = cor_test_wrapper, c("MIQ", "RAT"), size = size)
  cor_plot <- distribution_plot(boot_cor, 
                                var = "estimate", 
                                raw_value = raw_cor, 
                                title = "Pearson correlation coefficient MIQ ~ RAT")
  messagef("Mean correlation with %d bootstrap samples: %.3f (+/- %.3f), raw correlation: %.3f", size, mean(boot_cor$estimate, na.rm = T), sd(boot_cor$estimate, na.rm = T)/sqrt(size),raw_cor)
  print(cor_plot)
  return()
  invisible(readline(prompt="Press [enter] to continue"))
  #2. Simple linear regression
  #2.1 R^2
  messagef("Running linear regression bootstrap with %d iterations for adj. R^2", size)
  
  raw_lm <- data %>% lm_wrapper(broom_FUN = broom::glance)
  boot_lm <- bootstrapper(data, FUN = lm_wrapper, c("MIQ", "RAT"), size = size, broom_FUN = broom::glance)
  r_adj_plot <- distribution_plot(boot_lm, 
                                  var = "adj.r.squared", 
                                  raw_value = raw_lm$adj.r.squared, 
                                  title = "adj. R^2: academic.overall ~ MIQ + RAT")
  messagef("Mean Adj. r^2 with %d bootstrap samples: %.3f (+/- %.3f), raw adj. R^2: %.3f", size, mean(boot_lm$adj.r.squared, na.rm = T), sd(boot_lm$adj.r.squared, na.rm = T)/sqrt(size), raw_lm$adj.r.squared)
  print(r_adj_plot)
  invisible(readline(prompt="Press [enter] to continue"))
  
  
  #2.2 beta weights
  messagef("Running linear regression bootstrap with %d iterations for beta weights", size)
  raw_lm <- data %>% lm_wrapper(broom_FUN = broom::tidy)
  raw_beta_MIQ <- raw_lm %>% filter(term == "MIQ.score") %>% pull(estimate)
  raw_beta_RAT <- raw_lm %>% filter(term == "RAT.score") %>% pull(estimate)
  
  boot_lm <- bootstrapper(data, FUN = lm_wrapper, c("MIQ", "RAT"), size = size, broom_FUN = broom::tidy)
  beta_MIQ_plot <- distribution_plot(boot_lm %>% filter(term == "MIQ.score"), 
                                     var = "estimate", 
                                     raw_value = raw_beta_MIQ, 
                                     title = "beta_MIQ")
  messagef("Mean beta_MIQ with %d bootstrap samples: %.3f (+/- %.3f), raw beta_MIQ: %.3f", 
           size, 
           mean(boot_lm %>% filter(term == "MIQ.score") %>% pull(estimate), 
                na.rm = T), 
           sd(boot_lm %>% filter(term == "MIQ.score") %>% pull(estimate), na.rm = T)/sqrt(size), 
           raw_beta_MIQ)
  print(beta_MIQ_plot)
  invisible(readline(prompt="Press [enter] to continue"))
  
  beta_RAT_plot <- distribution_plot(boot_lm %>% filter(term == "RAT.score"), 
                                     var = "estimate", 
                                     raw_value = raw_beta_RAT, 
                                     title = "beta_RAT")
  messagef("Mean beta_RAT with %d bootstrap samples: %.3f (+/- %.3f), raw beta_RAT: %.3f", 
           size, 
           mean(boot_lm %>% filter(term == "RAT.score") %>% pull(estimate), 
                na.rm = T), 
           sd(boot_lm %>% filter(term == "RAT.score") %>% pull(estimate), na.rm = T)/sqrt(size), 
           raw_beta_RAT)  
  print(beta_RAT_plot)
  invisible(readline(prompt="Press [enter] to continue"))
  
  #3. simple SEM
  messagef("Running FA-SEM bootstrap with %d iterations", size)
  raw_sem <- data %>% sem_wrapper() 
  raw_betas <- raw_sem %>% 
    filter(lhs == "model_music", op == "=~") %>% 
    arrange(rhs) %>% pull(est.std)
  
  boot_sem<- bootstrapper(data, 
                          FUN = sem_wrapper, c("MDI", "MPT", "BAT"), 
                          size = size)
  sem_beta_plot <- distribution_plot(boot_sem %>% filter(lhs == "model_music", op == "=~"), 
                                     var = "est.std", 
                                     group_var = "rhs",
                                     raw_value = raw_betas, 
                                     title = "Factor loadings model_music ~ MDI + MPT + BAT")
  #browser()
  cmp_values <- bind_cols(
    boot_sem %>% 
      filter(lhs == "model_music", op == "=~") %>% 
      group_by(rhs) %>% 
      summarise(mean_value = mean(est.std, na.rm =  T)) %>% 
      ungroup(),
    tibble(
      raw_value = raw_betas)
  ) %>% rename(Predictor = rhs, Raw = raw_value, `Bootstrap Mean` = mean_value) 
  messagef("SEM betas with %d bootstrap samples and raw values", size)
  require(knitr)
  print(knitr::kable(cmp_values, digits = 3))
  print(sem_beta_plot)

  list(cor_plot = cor_plot, 
       r_adj_plot = r_adj_plot,
       beta_MIQ_plot = beta_MIQ_plot,
       beta_RAT_plot = beta_RAT_plot, 
       sem_beta_plot = sem_beta_plot) %>% invisible()
  
}
#demo_bootstrapper(size = 100)