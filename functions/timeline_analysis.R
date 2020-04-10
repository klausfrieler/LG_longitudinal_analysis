library(lme4)
library(lmerTest)
library(tidyverse)
library(broom)
library(broom.mixed)

#'Helper frunction that creates a version of data fraem  @data which each row repeated @n times
rep_rows_dfr <- function(data, n){
  map_dfr(1:nrow(data), function(i){
    map_dfr(1:n, function(x){
      data[i,]
    })
  })  
}

#' Helper function to bind data frames with different numbers of rows column-wise 
#' Needed to multiverse outputs
bind_cols_m <- function(tib1, tib2){
  if(nrow(tib2) < nrow(tib1)){
    tmp <- tib1
    tib1 <- tib2
    tib2 <- tmp
  }
  n1 <- nrow(tib1)
  n2 <- nrow(tib2)
  mult <- n2 %% n1
  if(mult != 0){
    stop("Rows of longer tibble not a multiple of shorter")
  }
  mult <- n2 %/% n1 
  tib1 <- rep_rows_dfr(tib1, mult)
  bind_cols(tib1, tib2)

}

#'Creates a filtered version of the original data for multiverse analysis
#'
compile_data <- function(
  raw_data, # input data goes here
  include_only_best_shot, 
  exclude_hearing_impaired,
  remove_2015,
  only_standard_gender
) {
  df <- raw_data
  if (include_only_best_shot) df <- df %>% filter(best_shot == T)
  if (exclude_hearing_impaired) df <- df %>% filter(hearing_impairment.binary == F)
  if (remove_2015) df <- df %>% filter(test_year != 2015)
  if (only_standard_gender) df <- df %>% filter(!is.na(gender), gender %in% c("Female", "Male"))
  df
}

#'The actual analysis, here a sinple linear mixed regression of music perception over age with random
#'intercepts per participant
#'
timelime_analysis <- function(data){
  mod1 <- lmer(music_perception ~ age + (1|p_id), data = data, REML = TRUE)
  #browser()
  list(
    coef = tidy(mod1),
    fit = glance(mod1),
    n_data = tibble(n_data = nrow(data))
  )
}
#'Applies a certain filter @compilation_data 
#'and runs whatever function @fun you like of the different data views of @raw_data
#'
run_function_with_compilation_parameters <- function(
  fun,
  raw_data, 
  compilation_parameters
) {
  stopifnot(is.data.frame(raw_data), is.list(compilation_parameters))
  args <- c(raw_data = list(raw_data), compilation_parameters)
  compiled_data = do.call(compile_data, args)
  fun(compiled_data)
}
#'Run funciotn @fun over all different combination of @possible_compilation_parameters
#'on @raw_data
#'
run_function_with_all_compilation_parameters <- function(
  fun,
  raw_data,
  possible_compilation_parameters
) {
  df <- do.call(expand_grid, possible_compilation_parameters)
  df$res <- pmap(df, function(...) run_function_with_compilation_parameters(
    fun,
    raw_data, 
    compilation_parameters = list(...)
  ))
  df
}
#'Helper funtion to add an ID for a specific multiverse parameter setting to @data, 
#'helpful for further analysis of multiverse results
#'
add_multi_id <- function(data){
  data %>% mutate(multi_id = sprintf("%d%d%d%d",
                                     include_only_best_shot, 
                                     exclude_hearing_impaired, 
                                     remove_2015, only_standard_gender))  
}
#'Appliation of timeline_analyis on the @data using a multiverse of
#'eight possible parameter combination (best_shot T/F, hearing_impairment.binary T/F,
#'with or withou 2015 and using standardized (non-p.c.) gender). 
#'
#'Note: This is only an demostration example!
#'
multiverse_output <- function(data){
  ret <- 
    run_function_with_all_compilation_parameters(
      timelime_analysis,
      data,
      list(
        include_only_best_shot = c(FALSE, TRUE),
        exclude_hearing_impaired = c(FALSE, TRUE),
        remove_2015 = c(FALSE, TRUE),
        only_standard_gender = c(FALSE, TRUE)
      ))
  #browser()
  
  coefs  <- bind_cols_m(ret %>% select(1:4), map_dfr(ret %>% pull(res), function(x) x$coef))
  fits   <- bind_cols_m(ret %>% select(1:4), map_dfr(ret %>% pull(res), function(x) x$fit))
  n_data <- bind_cols_m(ret %>% select(1:4), map_dfr(ret %>% pull(res), function(x) x$n_data))
  
  coefs  <- coefs %>% add_multi_id()
  fits   <- fits %>% add_multi_id()
  n_data <- n_data %>%  add_multi_id()
  coefs  <- coefs %>% left_join(n_data %>% select(multi_id, n_data), by = "multi_id")
  fits   <- fits %>% left_join(n_data %>% select(multi_id, n_data), by = "multi_id")
  list(coefs = coefs, fits = fits)
}