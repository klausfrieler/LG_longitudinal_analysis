library(lavaan)
library(ggplot2)
library(tidyverse)

#'variables of interest
FA_variables <- c("p_id", 
                  "country", 
                  "school", 
                  "gender", 
                  "age", 
                  "BAT.score", 
                  "MDI.score", 
                  "MPT.score", 
                  "best_shot", 
                  "hearing_impairment.binary")
                           

# Basic FA model
model_music <-
  '# measurement model
  model_music =~ MDI.score + MPT.score + BAT.score'

#' help function to add test year to variable names for aggregation
expand_variable_names <-  function(data, year){
  nd <- which(names(data) %in% setdiff(names(data), c("p_id", "country", "school", "gender", "test_year")))
  names(data)[nd] <- sprintf("%s.%s", names(data)[nd], year) 
  data
} 

#' Do some nicht factor analyses
do_factor_analysis <- function(data = master, year, variables = FA_variables ){
  data <- data %>% filter(test_year == year) %>% select(variables)
  browser()
  fit_model_music <- sem(model_music, data = data, missing = "fiml", estimator = "mlr")
  #summary <- standardizedsolution(fit_model_music, standardized = TRUE, fit.measures = TRUE,  )
  fscores <- lavPredict(fit_model_music, assemble = TRUE)
  #browser()
  # add factor scores and mean scores to the dataset
  data$fscores <- as.vector(fscores[, 1])
  data$meanscores <- rowMeans(data %>% select(BAT.score, MDI.score, MPT.score))
  data$test_year <- year
  return(data)
}

#' setup the laten growth analysis, 
#' adds ds_test_year and ds_age and long_format_fscores to global workspace
setup_latent_growth_analysis <- function(){
  #' select cases in the dataset where participants were  tested more than once
  good_ids <- master %>% 
    filter(test_year %in% c(2017, 2018, 2019)) %>% 
    group_by(p_id) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>% 
    pull(p_id) %>% 
    unique()
  
  # do and aggregrate factor analyses
  ds_2017 <- do_factor_analysis(master %>% filter(p_id %in% good_ids), 2017) 
  ds_2018 <- do_factor_analysis(master %>% filter(p_id %in% good_ids), 2018) 
  ds_2019 <- do_factor_analysis(master %>% filter(p_id %in% good_ids), 2019) 
  
  # merge data
  ds_test_year <- ds_2017 %>% expand_variable_names(2017)  %>% select(-test_year) %>%
    left_join(ds_2018 %>% expand_variable_names(2018) %>% select(-country, -school, -gender, -test_year), by = "p_id") %>% 
    left_join(ds_2019 %>% expand_variable_names(2019) %>% select(-country, -school, -gender, -test_year), by = "p_id")
  
  
  # end with 1195 cases
  
  # match data according to participant age
  
  long_format_fscores <- rbind(ds_2017, ds_2018, ds_2019)
  ds_age <- long_format_fscores %>% 
    pivot_wider(id_cols = c("p_id","country","school","gender","best_shot","hearing_impairment.binary"), 
                names_from = "age", 
                values_from = c("BAT.score", "MDI.score", "MPT.score", "fscores", "meanscores"))
  
  # delete cases in the dataset where participants were only tested once
  missingvalues <- ds_age %>% select(starts_with(c("BAT", "MDI", "MPT"))) %>% is.na() 
  keep_cases <- rowSums(missingvalues) <= 21
  # begin with 5133 cases
  ds_age <- ds_age %>% filter(keep_cases)
  # end with 4734 cases
  assign("ds_test_year", ds_test_year, globalenv())
  assign("ds_age", ds_age, globalenv())
  assign("long_format_fscores", long_format_fscores, globalenv())
}

#' Some nice plots
make_diagnostic_plots <- function(data = long_format_fscores){
  ggplot(data, aes(x = age, y = fscores)) + geom_boxplot(aes(group=age))
  ggsave('figs/boxplot-age-fscores.jpg')
  ggplot(data, aes(x = age, y = meanscores)) + geom_boxplot(aes(group=age))
  ggsave('figs/boxplot-age-meanscores.jpg')
  ggplot(data, aes(x = test_year, y = fscores)) + geom_boxplot(aes(group = test_year))
  ggsave('figs/boxplot-test-year-fscores.jpg')
  ggplot(data, aes(x = test_year, y = meanscores)) + geom_boxplot(aes(group = test_year))
  ggsave('figs/boxplot-test-year-meanscores.jpg') 
}
#' helper function for running all the laten growth models
#' retrieves and modifies datasets
get_dataset <- function(data_set){
  if(data_set  == "ds_test_year" || data_set == "ds_age"){
    return(get(data_set))
  }
  if(data_set == "rescaled_ds_test_year"){
    tmp <- ds_test_year
    tmp$fscores.2017 <- ((tmp$fscores.2017 + 2) / 6)*100
    tmp$fscores.2018 <- ((tmp$fscores.2018 + 2) / 6)*100
    tmp$fscores.2019 <- ((tmp$fscores.2019 + 2) / 6)*100
    return(tmp)
  }
  if(dataset == "adjusted_ds_test_year"){
    fit3_m1 <- 
      lavaan(models[["music_3_m1"]]$model, 
          data = models[["music_3_m1"]]$data_set, 
          orthogonal = T, 
          missing = "ML", 
          estimator = "MLR")
    fscores <- lavPredict(fit3_m1, assemble = TRUE)
    tmp <- ds_test_year
    tmp$fscores.2017 <- fscores [,1]
    tmp$fscores.2018 <- fscores [,2]
    tmp$fscores.2019 <- fscores [,3]
    return(tmp)
  }
}
#' Loads model definition from latent_growth_models.R
#' and runs all the models
#' returns list of model fits and stuff
#'  
run_latent_growth_models <- function(){
  source("./scripts/latent_growth_models.R")
  map(names(models), function(x){
    tryCatch({
      fit <- lavaan(models[[x]]$model, 
                  data = get_dataset(models[[x]]$data_set), 
                  orthogonal = T, 
                  missing = "ML", 
                  estimator= "MLR")
      },
      error = function(e){
        return(NULL)
      },
      warning = function(e){
        messagef("Warning for '%s'", x)
        warnings()
      })

    list(name = x,
         description = models[[x]]$description,
         fit = standardizedSolution(fit), 
         fit_measures = fitmeasures(fit, c("chisq", "df", "pvalue", "rmsea", "tli", "cfi")))
  }) %>% setNames(names(models)) %>% compact()
}



