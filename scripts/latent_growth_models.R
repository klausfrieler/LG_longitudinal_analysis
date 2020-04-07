## factorial invariance between test-year (longitudinal CFA) (Mair, 2018, pp. 52-55; Little, 2013, pp.137-179)

# 1. specify the longitudinal mull model (Little, 2013, pp.112-118 & 160)
models <- list()

models[["music_0"]] <- list(
  descripton = "specify the longitudinal mull model (Little, 2013, pp.112-118 & 160)",
  model = '
   BAT.score.2017 ~~ V1*BAT.score.2017
   BAT.score.2018 ~~ V1*BAT.score.2018
   BAT.score.2019 ~~ V1*BAT.score.2019
   MDI.score.2017 ~~ V2*MDI.score.2017
   MDI.score.2018 ~~ V2*MDI.score.2018
   MDI.score.2019 ~~ V2*MDI.score.2019
   MPT.score.2017 ~~ V3*MPT.score.2017
   MPT.score.2018 ~~ V3*MPT.score.2018
   MPT.score.2019 ~~ V3*MPT.score.2019
   
   BAT.score.2017 ~ T1*1
   MDI.score.2017 ~ T2*1
   MPT.score.2017 ~ T3*1
   BAT.score.2018 ~ T1*1
   MDI.score.2018 ~ T2*1
   MPT.score.2018 ~ T3*1
   BAT.score.2019 ~ T1*1
   MDI.score.2019 ~ T2*1
   MPT.score.2019 ~ T3*1
  ', 
  data_set = "ds_test_year")


# 2. test for configural invariance
models[["music_1"]] <- list(
  description = "test for configural invariance",
  model = '
   # defining factors
   music.2017 =~ L1*BAT.score.2017 + MDI.score.2017 + MPT.score.2017
   music.2018 =~ BAT.score.2018 + MDI.score.2018 + MPT.score.2018
   music.2019 =~ BAT.score.2019 + MDI.score.2019 + MPT.score.2019
  
   # error variances
   BAT.score.2017 ~~ BAT.score.2017
   BAT.score.2018 ~~ BAT.score.2018
   BAT.score.2019 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2017
   MDI.score.2018 ~~ MDI.score.2018
   MDI.score.2019 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2017
   MPT.score.2018 ~~ MPT.score.2018
   MPT.score.2019 ~~ MPT.score.2019
   
   # latent factor variances
   music.2017 ~~ music.2017
   music.2018 ~~ music.2018
   music.2019 ~~ music.2019
   
   # covariances of errors
   BAT.score.2017 ~~ BAT.score.2018 + BAT.score.2019
   BAT.score.2018 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2018 + MDI.score.2019
   MDI.score.2018 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2018 + MPT.score.2019
   MPT.score.2018 ~~ MPT.score.2019
   
   # covariances of latent facors
   music.2017 ~~ music.2018 + music.2019
   music.2018 ~~ music.2019
   
   # latent means (not estimated)
   music.2017 ~ NA*1
   music.2018 ~ NA*1
   music.2019 ~ NA*1
   
   ## intercepts
   BAT.score.2017 ~ T1*1
   BAT.score.2018 ~ T4*1
   BAT.score.2019 ~ T7*1
   MDI.score.2017 ~ T2*1
   MDI.score.2018 ~ T5*1
   MDI.score.2019 ~ T8*1
   MPT.score.2017 ~ T3*1
   MPT.score.2018 ~ T6*1
   MPT.score.2019 ~ T9*1
  ',
  data_set = "ds_test_year")



# 3. test for weak invariance
models[["music_2"]] <-list(
  description = "test for weak invariance",
  model = '
   # defining factors
   music.2017 =~ L1*BAT.score.2017 + L2*MDI.score.2017 + L3*MPT.score.2017
   music.2018 =~ L1*BAT.score.2018 + L2*MDI.score.2018 + L3*MPT.score.2018
   music.2019 =~ L1*BAT.score.2019 + L2*MDI.score.2019 + L3*MPT.score.2019
   
   # error variances
   BAT.score.2017 ~~ BAT.score.2017
   BAT.score.2018 ~~ BAT.score.2018
   BAT.score.2019 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2017
   MDI.score.2018 ~~ MDI.score.2018
   MDI.score.2019 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2017
   MPT.score.2018 ~~ MPT.score.2018
   MPT.score.2019 ~~ MPT.score.2019
   
   # latent factor variances
   music.2017 ~~ music.2017
   music.2018 ~~ music.2018
   music.2019 ~~ music.2019
   
   # covariances of errors
   BAT.score.2017 ~~ BAT.score.2018 + BAT.score.2019
   BAT.score.2018 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2018 + MDI.score.2019
   MDI.score.2018 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2018 + MPT.score.2019
   MPT.score.2018 ~~ MPT.score.2019
   
   # covariances of latent facors
   music.2017 ~~ music.2018 + music.2019
   music.2018 ~~ music.2019
   
   # latent means (not estimated)
   music.2017 ~ 1
   music.2018 ~ 1
   music.2019 ~ 1
   
   ## intercepts
   BAT.score.2017 ~ T1*1
   BAT.score.2018 ~ T4*1
   BAT.score.2019 ~ T7*1
   MDI.score.2017 ~ T2*1
   MDI.score.2018 ~ T5*1
   MDI.score.2019 ~ T8*1
   MPT.score.2017 ~ T3*1
   MPT.score.2018 ~ T6*1
   MPT.score.2019 ~ T9*1
  ',
  data_set = "ds_test_year")



# 4. test for strong invariance

models[["music_3"]] <- list(
  description = "test for strong invariance",
  model = '
   # defining factors
   music.2017 =~ L1*BAT.score.2017 + L2*MDI.score.2017 + L3*MPT.score.2017
   music.2018 =~ L1*BAT.score.2018 + L2*MDI.score.2018 + L3*MPT.score.2018
   music.2019 =~ L1*BAT.score.2019 + L2*MDI.score.2019 + L3*MPT.score.2019
   
   # error variances
   BAT.score.2017 ~~ BAT.score.2017
   BAT.score.2018 ~~ BAT.score.2018
   BAT.score.2019 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2017
   MDI.score.2018 ~~ MDI.score.2018
   MDI.score.2019 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2017
   MPT.score.2018 ~~ MPT.score.2018
   MPT.score.2019 ~~ MPT.score.2019
   
   # latent factor variances
   music.2017 ~~ music.2017
   music.2018 ~~ music.2018
   music.2019 ~~ music.2019
   
   # covariances of errors
   BAT.score.2017 ~~ BAT.score.2018 + BAT.score.2019
   BAT.score.2018 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2018 + MDI.score.2019
   MDI.score.2018 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2018 + MPT.score.2019
   MPT.score.2018 ~~ MPT.score.2019
   
   # covariances of latent facors
   music.2017 ~~ music.2018 + music.2019
   music.2018 ~~ music.2019
   
   # latent means (not estimated)
   music.2017 ~ 1
   music.2018 ~ 1
   music.2019 ~ 1
   
   ## intercepts
   BAT.score.2017 ~ T1*1
   BAT.score.2018 ~ T1*1
   BAT.score.2019 ~ T1*1
   MDI.score.2017 ~ T2*1
   MDI.score.2018 ~ T2*1
   MDI.score.2019 ~ T2*1
   MPT.score.2017 ~ T3*1
   MPT.score.2018 ~ T3*1
   MPT.score.2019 ~ T3*1
  ',
  data_set = "ds_test_year")


## test for strong invariance didn't pass

# 5. partial invariance
## 5-1. take a look at the modification indices
#mod_music3 <- modindices(fit3, free.remove = FALSE)
#music3_intercepts <- mod_music3 [mod_music3$op == "~1",]
#na.omit(music3_intercepts [order(mod_music3$mi, decreasing = TRUE), ])
## MPT.score.2017; MDI.score.2017; MDI.score.2019; MPT.score.2019 produce largest misfit

models[["music_3_m1"]] <- list(
  model = '
   # defining factors
   music.2017 =~ L1*BAT.score.2017 + L2*MDI.score.2017 + L3*MPT.score.2017
   music.2018 =~ L1*BAT.score.2018 + L2*MDI.score.2018 + L3*MPT.score.2018
   music.2019 =~ L1*BAT.score.2019 + L2*MDI.score.2019 + L3*MPT.score.2019
   
   # error variances
   BAT.score.2017 ~~ BAT.score.2017
   BAT.score.2018 ~~ BAT.score.2018
   BAT.score.2019 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2017
   MDI.score.2018 ~~ MDI.score.2018
   MDI.score.2019 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2017
   MPT.score.2018 ~~ MPT.score.2018
   MPT.score.2019 ~~ MPT.score.2019
   
   # latent factor variances
   music.2017 ~~ music.2017
   music.2018 ~~ music.2018
   music.2019 ~~ music.2019
   
   # covariances of errors
   BAT.score.2017 ~~ BAT.score.2018 + BAT.score.2019
   BAT.score.2018 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2018 + MDI.score.2019
   MDI.score.2018 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2018 + MPT.score.2019
   MPT.score.2018 ~~ MPT.score.2019
   
   # covariances of latent facors
   music.2017 ~~ music.2018 + music.2019
   music.2018 ~~ music.2019
   
   # latent means (not estimated)
   music.2017 ~ 1
   music.2018 ~ 1
   music.2019 ~ 1
   
   ## intercepts
   BAT.score.2017 ~ T1*1
   BAT.score.2018 ~ T1*1
   BAT.score.2019 ~ T1*1
   MDI.score.2017 ~ T2*1
   MDI.score.2018 ~ T2*1
   MDI.score.2019 ~ T2*1
   MPT.score.2017 ~ T3*1
   MPT.score.2018 ~ T4*1
   MPT.score.2019 ~ T4*1
  ',
  data_set = "ds_test_year")



# 6. extract factor scores from this partial invariance model


## latent growth model (using the partial measurement model above)
models[["growth_model"]] <- list(
  description = "extract factor scores from this partial invariance model",
  model = '
   ## measurement model
   # defining factors
   music.2017 =~ L1*BAT.score.2017 + L2*MDI.score.2017 + L3*MPT.score.2017
   music.2018 =~ L1*BAT.score.2018 + L2*MDI.score.2018 + L3*MPT.score.2018
   music.2019 =~ L1*BAT.score.2019 + L2*MDI.score.2019 + L3*MPT.score.2019
   
   # error variances
   BAT.score.2017 ~~ BAT.score.2017
   BAT.score.2018 ~~ BAT.score.2018
   BAT.score.2019 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2017
   MDI.score.2018 ~~ MDI.score.2018
   MDI.score.2019 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2017
   MPT.score.2018 ~~ MPT.score.2018
   MPT.score.2019 ~~ MPT.score.2019
   
   # latent factor variances
   music.2017 ~~ music.2017
   music.2018 ~~ music.2018
   music.2019 ~~ music.2019
   
   # covariances of errors
   BAT.score.2017 ~~ BAT.score.2018 + BAT.score.2019
   BAT.score.2018 ~~ BAT.score.2019
   MDI.score.2017 ~~ MDI.score.2018 + MDI.score.2019
   MDI.score.2018 ~~ MDI.score.2019
   MPT.score.2017 ~~ MPT.score.2018 + MPT.score.2019
   MPT.score.2018 ~~ MPT.score.2019
   
   # covariances of latent facors
   music.2017 ~~ music.2018 + music.2019
   music.2018 ~~ music.2019
   
   # latent means (not estimated)
   music.2017 ~ 1
   music.2018 ~ 1
   music.2019 ~ 1
   
   # intercepts
   BAT.score.2017 ~ T1*1
   BAT.score.2018 ~ T1*1
   BAT.score.2019 ~ T1*1
   MDI.score.2017 ~ T2*1
   MDI.score.2018 ~ T2*1
   MDI.score.2019 ~ T2*1
   MPT.score.2017 ~ T3*1
   MPT.score.2018 ~ T4*1
   MPT.score.2019 ~ T4*1
   
   ## intercept and slope with fixed coefficients
   i =~ 1*music.2017 + 1*music.2018 + 1*music.2019
   s =~ 0*music.2017 + 1*music.2018 + 2*music.2019
  ',
  data_set = "adjusted_ds_test_year")

#fscores <- lavPredict(fit3m1, assemble = TRUE)
#ds_test_year$fscores.2017 <- fscores [,1]
#ds_test_year$fscores.2018 <- fscores [,2]
#ds_test_year$fscores.2019 <- fscores [,3]


# model was not estimated correctly (slopes weren't estimated)

## latent growth model (using fscores extracted from the partial invariant model)
# rescale fscores
#ds_test_year$fscores.2017 <- ((ds_test_year$fscores.2017 + 2) / 6)*100
#ds_test_year$fscores.2018 <- ((ds_test_year$fscores.2018 + 2) / 6)*100
#ds_test_year$fscores.2019 <- ((ds_test_year$fscores.2019 + 2) / 6)*100

models[["model_fscores"]] <- list(
  descrition = "latent growth model (using fscores extracted from the partial invariant model)",
  model = '
    # intercept and slope with fixed coefficients
  i =~ 1*fscores.2017 + 1*fscores.2018 + 1*fscores.2019
  s =~ 0*fscores.2017 + 1*fscores.2018 + 2*fscores.2019',
  data_set = "rescaled_ds_test_year")


models[["model_meanscores"]] <- list(
  model = '
    # intercept and slope with fixed coefficients
    i =~ 1*meanscores.2017 + 1*meanscores.2018 + 1*meanscores.2019
    s =~ 0*meanscores.2017 + 1*meanscores.2018 + 2*meanscores.2019',
  data_set = "rescaled_ds_test_year")


#latent growth model (ds_age)
models[["model_age"]] <- list(
  description = "latent growth model (over age)",
  model = '
    # intercept and slope with fixed coeff icients
    i =~ 1*fscores.10 + 1*fscores.11 + 1*fscores.12 + 1*fscores.13 + 1*fscores.14 + 1*fscores.15 + 1*fscores.16 + 1*fscores.17
    s =~ 0*fscores.10 + 1*fscores.11 + 2*fscores.12 + 3*fscores.13 + 4*fscores.14 + 5*fscores.15 + 6*fscores.16 + 7*fscores.17',
  data_set = "ds_age")



#save.image("2020-03-27-paperthon_working_progress.RData")

