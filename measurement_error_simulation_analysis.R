library(tidyverse)
library(log4r)
library(tictoc)
logger <- create.logger()
logfile(logger) <- "me_simu.log"
level(logger) <- "INFO"

version <- "0.4.0"
aws3_dir_simulations <- "measurement_error_study/simulations"
aws3_bucket <- "longgold.gold-msi.org"
bar <- paste(rep("-", 60), collapse ="")

#This hack avoids including library(mclust), due to bug in mclust
mclustBIC <- mclust::mclustBIC

source("setup.R")


messagef <- function(...){
  message(sprintf(...))
  log4r::info(logger, sprintf(...))
}

deco_messagef <- function(...){
  messagef(bar)
  messagef(...)
  messagef(bar)
}

download_simu_summary_data <- function(version = "0.4.0"){
  browser()
  aws.set_credentials("kf_aws-credentials.txt")
  summary_stats_files <- aws.s3::get_bucket_df(aws3_bucket, prefix = aws3_dir_simulations) %>%
    filter(str_detect(Key, version))  %>% 
    pull(Key)
  stats <- map_dfr(summary_stats_files, function(x){
    print(gsub("=", "%3D", x))
    tmp <- aws.s3::s3readRDS(object = x, bucket = aws3_bucket) 
    tmp$stats %>% mutate(error_level_f = factor(error_level, levels = c("none", "small", "medium", "large")))
  })
  pool <- map_dfr(summary_stats_files, function(x){
    print(gsub("=", "%3D", x))
    tmp <- aws.s3::s3readRDS(object = x, bucket = aws3_bucket) 
    tmp$pool %>% mutate(error_level_f = factor(error_level, levels = c("none", "small", "medium", "large")))
  })
  list(pool = pool, stats = stats)  
}

read_simu_summary_data <- function(dir = "simulations", version = "0.4.0"){
  files <- list.files(dir, pattern = "rds", full.names = T)
  browser()
  summary_stats_files <- files[str_detect(files, version)]
  if(length(summary_stats_files) == 0){
    stop("No files found")
  }
  stats <- map_dfr(summary_stats_files, function(x){
    tmp <- readRDS(x)
    browser()
    tmp$stats %>% mutate(error_level_f = factor(error_level, levels = c("none", "small", "medium", "large")))
  })
  pool <- map_dfr(summary_stats_files, function(x){
    tmp <- readRDS(x)
    tmp$pool %>% mutate(error_level_f = factor(error_level, levels = c("none", "small", "medium", "large")))
  })
  list(pool = pool, stats = stats)  
}

plot_summary_stats <- function(sum_stats, variable = "sum_stats", metrics = "rel_diff", with_na = T, yline = 0){
  sum_stats <- sum_stats %>% 
    filter(term == variable, with_na == !!with_na, type != "true") %>% 
    group_by(type, error_level_f, method, size) %>% 
    summarise(!!sym(metrics) := mean(!!sym(metrics)), se := se(!!sym(metrics)), .groups = "drop")

  q <- sum_stats %>% ggplot(aes(x = type, y = !!sym(metrics), colour = error_level_f)) 
  q <- q + geom_point() 
  q <- q + facet_grid(size ~ method) 
  q <- q + geom_hline(yintercept = yline) 
  q <- q + theme_bw() 
  q <- q + theme(legend.position = "top")
  q <- q + geom_line(aes(group = error_level_f)) 
  if(metrics %in% c("rel_diff", "abs_rel_diff")){
    q <- q + scale_y_continuous(labels = scales::percent) 
  }
  q <- q + ggtitle(sprintf("%s (NA = %s)", variable, with_na)) 
  q <- q + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  q <- q + geom_errorbar(aes(ymin = !!sym(metrics) - se, ymax = !!sym(metrics) + se), width = 1)
  q <- q + labs(colour = "Error Level")
  q
}

plot_all_summary_stats <-function(sum_stats, out_dir = "figs/simulations", with_error_level_none = T, no_na = F){
  metrics <- c( "rel_diff", "abs_rel_diff", "abs_diff", "attenuation")
  ylines <- c(0, 0, 0, 1)
  variable <- c("sum_stats", "MIQ.score", "MT.score", "MHE.score")
  comb <- tibble(metric = metrics, yline = ylines)
  comb <- cross_add_col(comb, variable)
  if(no_na){
    comb  <- cross_add_col(comb, with_na)
    
  }
  else{
    comb$with_na  <- TRUE
  }
  if(!with_error_level_none){
    sum_stats <- filter(errole_level != "none")  
  }
  for(r in 1:nrow(comb)){
    file_name <- file.path(out_dir, sprintf("%s-%s-NA=%s.png", comb[r, ]$variable, comb[r, ]$metric, comb[r, ]$with_na ))
    q <- plot_summary_stats(sum_stats, comb[r, ]$variable, comb[r, ]$metric, comb[r, ]$with_na, comb[r, ]$yline)
    ggsave(filename = file_name, plot = q)
  }
}