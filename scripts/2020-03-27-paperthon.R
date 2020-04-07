setwd("~/10_Paperthon_Mar2020/data/")
paperthon_data_wide <- read.csv("~/10_Paperthon_Mar2020/data/paperthon_data_wide.csv", sep=";")
paperthon_data_long <- read.csv("~/10_Paperthon_Mar2020/data/paperthon_data_long.csv", sep=";")
long_format <- paperthon_data_long [which(is.na(paperthon_data_long$test_year) %in% FALSE),]
selected_variables <- c("p_id", "country", "school", "gender", "age", "BAT.score", "MDI.score", "MPT.score", "best_shot", "hearing_impairment.binary")
factor_analysis_2017 <- long_format [long_format$test_year == 2017, selected_variables]
factor_analysis_2018 <- long_format [long_format$test_year == 2018, selected_variables]
factor_analysis_2019 <- long_format [long_format$test_year == 2019, selected_variables]                           

save.image("2020-03-27-paperthon_working_progress.RData")

library(lavaan)

# calculate factor scores 2017
model_music <-
'# measurement model
model_music =~ MDI.score + MPT.score + BAT.score'
fit_model_music <- sem(model_music, data=factor_analysis_2017, missing="fiml", estimator="mlr")
summary_2017 <- summary(fit_model_music, standardized=TRUE, fit.measures=TRUE)
fscores_2017 <- lavPredict(fit_model_music, assemble = TRUE)

# add factor scores and mean scores to the dataset 2017
factor_analysis_2017$fscores <- NULL
factor_analysis_2017$fscores <- fscores_2017
factor_analysis_2017$meanscroes <- NULL
factor_analysis_2017$meanscores <- rowMeans(factor_analysis_2017[,which(selected_variables %in% c("BAT.score","MDI.score","MPT.score"))])

# calculate factor scores 2018
fit_model_music <- sem(model_music, data=factor_analysis_2018, missing="fiml", estimator="mlr")
summary_2018 <- summary(fit_model_music, standardized=TRUE, fit.measures=TRUE)
fscores_2018 <- lavPredict(fit_model_music, assemble = TRUE)

# add factor scores and mean scores to the dataset 2018
factor_analysis_2018$fscores <- NULL
factor_analysis_2018$fscores <- fscores_2018
factor_analysis_2018$meanscroes <- NULL
factor_analysis_2018$meanscores <- rowMeans(factor_analysis_2018[,which(selected_variables %in% c("BAT.score","MDI.score","MPT.score"))])

# calculate factor scores 2019
fit_model_music <- sem(model_music, data=factor_analysis_2019, missing="fiml", estimator="mlr")
summary_2019 <- summary(fit_model_music, standardized=TRUE, fit.measures=TRUE)
fscores_2019 <- lavPredict(fit_model_music, assemble = TRUE)

# add factor scores and mean scores to the dataset 2019
factor_analysis_2019$fscores <- NULL
factor_analysis_2019$fscores <- fscores_2019
factor_analysis_2019$meanscroes <- NULL
factor_analysis_2019$meanscores <- rowMeans(factor_analysis_2019[,which(selected_variables %in% c("BAT.score","MDI.score","MPT.score"))])

# match data according to test year
selected_variables2 <- c("p_id","fscores","meanscores","best_shot", "hearing_impairment.binary","age","country","gender","school")
ds_2017 <- factor_analysis_2017[, which(names(factor_analysis_2017)%in% selected_variables2)]
ds_2018 <- factor_analysis_2018[, which(names(factor_analysis_2018)%in% selected_variables2)]
ds_2019 <- factor_analysis_2019[, which(names(factor_analysis_2019)%in% selected_variables2)]
names(ds_2017)[2:length(names(ds_2017))] <- paste(names(ds_2017)[2:length(names(ds_2017))], "2017", sep = ".")
names(ds_2018)[2:length(names(ds_2018))] <- paste(names(ds_2018)[2:length(names(ds_2018))], "2018", sep = ".")
names(ds_2019)[2:length(names(ds_2019))] <- paste(names(ds_2019)[2:length(names(ds_2019))], "2019", sep = ".")
# merge data
ds_1718 <- merge(ds_2017, ds_2018, by = "p_id", all = T)
ds_test_year <- merge(ds_1718, ds_2019, by = "p_id", all = T)

# match data according to participant age
ds_2017 <- factor_analysis_2017[, which(names(factor_analysis_2017)%in% selected_variables2)]
ds_2018 <- factor_analysis_2018[, which(names(factor_analysis_2018)%in% selected_variables2)]
ds_2019 <- factor_analysis_2019[, which(names(factor_analysis_2019)%in% selected_variables2)]
ds_2017$test_year <- 2017
ds_2018$test_year <- 2018
ds_2019$test_year <- 2019
long_format_fscores <- rbind(ds_2017, ds_2018, ds_2019)
ds_age <- reshape(long_format_fscores, idvar = c("p_id","country","school","gender","best_shot","hearing_impairment.binary","test_year"),
                  timevar = "age", direction = "wide")

# plotting
library(ggplot2)

ggplot(long_format_fscores, aes(x = age, y = fscores)) + geom_boxplot(aes(group=age))
ggsave('boxplot-age-fscores.jpg')
ggplot(long_format_fscores, aes(x = age, y = meanscores)) + geom_boxplot(aes(group=age))
ggsave('boxplot-age-meanscores.jpg')
ggplot(long_format_fscores, aes(x = test_year, y = fscores)) + geom_boxplot(aes(group = test_year))
ggsave('boxplot-test-year-fscores.jpg')
ggplot(long_format_fscores, aes(x = test_year, y = meanscores)) + geom_boxplot(aes(group = test_year))
ggsave('boxplot-test-year-meanscores.jpg')

# latent grwoth model

model <- '
  # intercept and slope with fixed coefficients
i =~ 1*fscores.2017 + 1*fscores.2018 + 1*fscores.2019
s =~ 0*fscores.2017 + 1*fscores.2018 + 2*fscores.2019'

fit <- growth(model, data = ds_test_year, missing ="FIML", estimator="MLR")
summary(fit, standardized=TRUE)

model <- '
  # intercept and slope with fixed coefficients
i =~ 1*meanscores.2017 + 1*meanscores.2018 + 1*meanscores.2019
s =~ 0*meanscores.2017 + 1*meanscores.2018 + 2*meanscores.2019'

fit <- growth(model, data = ds_test_year, missing ="FIML", estimator="MLR")
summary(fit, standardized=TRUE)

