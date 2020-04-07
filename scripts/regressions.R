#### LongGold Paperthon Draft Analysis ###

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(lavaan)
library(blavaan)
library(tidyverse)
library(riclpmr)

#devtools::install_github('jflournoy/riclpmr')

###Model
run_lmer_model <- function(data = master){
  ### Model 1: Intraindividual change in musical ability
  ### Note that 'music_perception' is just the average across all test scores that an individual did in a particular year
  fit <- lmer(music_perception ~ age+(1|p_id), data = data, REML = TRUE)
  list(glance = broom.mixed::glance(fit), tidy = broom.mixed::tidy(fit))
} 

###Graph
lmer_plot <- function(data = master){
  g1 <- data %>% ggplot(aes(x = age, y = music_perception, group = p_id)) + geom_line(colour='grey') 
  g1 <- g1 + stat_summary(fun = mean, aes(group = 1), geom = "line", lwd = 3) 
  g1 <- g1 + scale_x_continuous(breaks = 10:17, name = "Age")
  g1 <- g1  + scale_y_continuous(name = "Average Listening test Score")
  g1
  
}

###Alternative Model using Bayesian approach and brms package
run_brms_regression <- function(data = master){
  mod1_brms <- brm(bf(MDI.score  ~ age + ( 1 | p_id)),
                   data = data, 
                   family = gaussian(),
                   prior = prior(normal(0, 1)),
                   cores = 4)
  mod1_brms_me <- brm(bf(MDI.score | se(MDI.error, sigma = TRUE) ~ age+(1|p_id)),
                      data = data, 
                      family = gaussian(),
                      prior = prior(normal(0, 1)),
                      cores = 4)
  
  list(brms = mod1_brms, brms_me = mod1_brms_me)
}

### Q4 ###


run_cross_lagged_sem <- function(data = master_wide){
  var_groups <- list(
    x = c("MIQ.score_11",  "MIQ.score_12",  "MIQ.score_13"),
    y = c("music_perception_11", "music_perception_12", "music_perception_13"))
  model_text <- riclpmr::riclpm_text(var_groups)
  #cat(model_text)
  riclpm_fit <- riclpmr::lavriclpm(riclpmModel = model_text, data = data, blavaan = T)
 
  blavaan_fit <- lavaan(model_text, data = data,
                missing = 'fiml', #for the missing data!
                int.ov.free = F,
                int.lv.free = F,
                auto.fix.first = F,
                auto.fix.single = F,
                auto.cov.lv.x = F,
                auto.cov.y = F,
                auto.var = F)
  list(riclpm_fit = standardizedSolution(riclpm_fit), 
       blavaan_fit = standardizedSolution(blavaan_fit))
}




