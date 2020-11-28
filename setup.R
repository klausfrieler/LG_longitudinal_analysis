library(tidyverse)
library(aws.s3)

source("functions/aws.R")
source("functions/data_wrangling.R")
source("functions/utils.R")

aws.set_credentials_if_unset()

#' Setups the workspace with all the data needed, pulling from S3.
#' Call this function always first to keep your data up-to-date. 
#' 
setup_workspace <- function(){

  messagef("Loading workspace...")
  
  data <- load_data()
  prepare_data(data) 
  #browser()
  assign("data_definitions", get_data_definitions(), globalenv())
  assign("raw_data", data, globalenv()) 
}


#' Prepares the data for this specific analysis project.
#' Filters for age, mainly
#' Overrides gender with a smoothed version of gender 
#' Prepares wide version of data with respect to age
#' Put objects "master" and "master_wide" in the global workspace for convenient use.
#' 
prepare_data <- function(data){
  master <- data %>% filter(age > 9, age < 18, p_id > 0, !is.na(age), !is.na(p_id))
  double_agers <- 
    master %>% 
    group_by(p_id) %>% 
    summarise(l = n_distinct(age), 
              n = length(age), 
              m = l == n) %>% 
    filter(!m) %>% 
    pull(p_id) 
  master <- master %>% filter(!(p_id %in% double_agers))
  master <- master %>% mutate(gender = patch_gender) %>% select(-patch_gender)
  
  master_wide <- master %>% 
    pivot_wider(
      id_cols = c(p_id, country, school, gender), 
      names_from = age,
      values_from = c(academic.overall, 
                      BAT.error, 
                      BAT.score, 
                      best_shot, 
                      CCM, 
                      CCM.weekly_num_school_music_lessons, 
                      CCM.weekly_after_school_music_lessons, 
                      EDT.score,
                      GMS.active_engagement,
                      GMS.emotions, 
                      GMS.general, 
                      GMS.musical_training, 
                      GMS.perceptual_abilities, 
                      GMS.singing_abilities,
                      hearing_impairment.binary, 
                      hearing_impairment.full_text,
                      MDI.error,
                      MDI.score, 
                      MHE.general_score, 
                      MIQ.error, 
                      MIQ.score, 
                      MPT.error, 
                      MPT.score,
                      PIAT.score, 
                      PIAT.error, 
                      RAT.score, 
                      RAT.error, 
                      SES.educational_degree,
                      music_perception))
  #browser()
  master_wide <- master_wide %>% 
    left_join(master %>% distinct(p_id, measurements), by ="p_id")
  
  master <- master %>% 
    select(p_id, 
           country, 
           school, 
           gender,
           test_year, 
           test_wave,
           age, 
           measurements,
           academic.overall, 
           BAT.error, 
           BAT.score, 
           best_shot, 
           CCM, 
           CCM.weekly_num_school_music_lessons, 
           CCM.weekly_after_school_music_lessons, 
           EDT.score,
           GMS.active_engagement,
           GMS.emotions, 
           GMS.general, 
           GMS.musical_training, 
           GMS.perceptual_abilities, 
           GMS.singing_abilities,
           hearing_impairment.binary, 
           hearing_impairment.full_text,
           MDI.error,
           MDI.score, 
           MHE.general_score, 
           MIQ.error, 
           MIQ.score, 
           MPT.error, 
           MPT.score,
           PIAT.score, 
           PIAT.error, 
           RAT.score, 
           RAT.error,
           SES.educational_degree,
           music_perception)
  assign("master", master, globalenv())
  assign("master_wide", master_wide, globalenv())

  messagef("Prepared and stored analysis data into 'master' (%d lines, %d distinct ids) 
          and 'master_wide' (%d lines, %d distinct ids)", nrow(master), length(unique(master$p_id)),
          nrow(master_wide), length(unique(master_wide$p_id)))
  if(n_distinct(master_wide$p_id) != nrow(master_wide)){
    messagef("Wide format malformed!")
    browser()
    stop()
  }  
}
#' Load the compiled data ("all.csv") from S3, 
#' adds some filtering and adds some standard vars and smaller tweaks
#' 
load_data <- function() {
  message("Loading data from S3...")
  data <- 
    aws.s3::s3read_using(
      function(x) 
        readr::read_csv(x, col_types = get_col_types()), 
      object = "data/compiled/all.csv",
      bucket = "longgold.gold-msi.org",
      opts = list(region = "eu-west-1"))
  #browser()
  messagef("...read %d lines of data\n", nrow(data))
  data[is.na(data$status) & !is.na(data$age),]$status <- "full"
  data[is.na(data$p_id), ]$p_id <- data[is.na(data$p_id), ]$p_id.old


  data$year_group <- make_year_group(data$year_group_uk, data$year_group_de, data$age, data$country)
  data$age_group <- make_age_group(data$year_group_uk, data$year_group_de, data$age, data$country)
  data <- data %>% group_by(p_id) %>% mutate(measurements = n()) %>% ungroup()
  data <- rescale_GMSI(data)
  data <- get_combined_var(data, c("BAT.score", "MPT.score", "MDI.score", "RAT.score"), 
                           label = "music_perception", scale  = F, is_error_var = F)  
  data <- get_combined_var(data, c("CCM", "DAC", "PAC") , label = "LAC")           
  data <- data %>% mutate(SES.educational_degree = 8 - SES.educational_degree)
  #data <- add_levelled_vars(data, aggregate = T)
  data
}

#' Read data definitions (variable names, types, use etc), needed for laoading data
#' 
get_data_definitions <- function(){
  read.csv("data/longgold_data_definition.csv", header = T, sep =";", stringsAsFactors = F) %>% 
    as_tibble()
}

#' Prepares column data types for correct loading of data
#' 
get_col_types <- function() {
  df <- get_data_definitions() %>% filter(!(variable %in% c("measurements", "year_group")))

  l <- df %>% 
    mutate(type = substr(type, 1, 1)) %>% 
    mutate(type = case_when(type == "f" ~ "c", 
                            type == "d" ~ "T", 
                            TRUE ~ type)) %>% 
    pull(type) %>% 
    as.list()
  names(l) <- df$variable
  do.call(readr::cols, l)
}

get_cross_section_version <- function(data){
  data %>% 
    group_by(p_id) %>% 
    mutate(latest = test_wave == max(test_wave)) %>% 
    ungroup() %>% 
    filter(latest)   
}

