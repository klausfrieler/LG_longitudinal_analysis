# Peter Harrison

source("./setup.R")
library(mgcv)
library(ggpubr)
library(scam)
theme_set(theme_pubr())

get_roughness <- function(par) {
  beta <- matrix(par, ncol = 1)
  # fda::bsplineS(breaks = BREAKS)
  basis <- fda::create.bspline.basis(breaks = BREAKS, norder = NORDER)
  penalty_matrix <- fda::bsplinepen(basis)
  roughness <- as.numeric(t(beta) %*% penalty_matrix %*% beta)
  roughness
}

get_response <- function(x, par) {
  beta <- matrix(par, ncol = 1)
  as.numeric(fda::bsplineS(x, breaks = BREAKS, norder = NORDER) %*% beta)
}

norm_par <- function(par) {
  anchors <- get_response(c(0, 1), par)
  scale <- anchors[2] - anchors[1]
  par / scale
}

# get_scaled_response <- function(x, par) {
#   anchors <- get_raw_response(c(0, 1), par)
#   shift <- anchors[1]
#   scale <- anchors[2] - anchors[1]
#   raw_response <- get_raw_response(x, par)
#   scaled_response <- (raw_response - shift) / scale
#   scaled_response
# }

plot_par <- function(par) {
  tibble(
    x = seq(from = -4, to = 4, length.out = 100),
    y = get_response(x, par)
  ) %>% 
    ggplot(aes(x, y)) + 
    geom_point() + 
    scale_x_continuous("Old IRT score") + 
    scale_y_continuous("New score") + 
    theme(aspect.ratio = 1)
}

rescale_error <- function(raw_score, raw_sem, par) {
  N <- 1e4
  raw_quartiles <- qnorm(c(0.25, 0.75), mean = raw_score, sd = raw_sem)
  new_quartiles <- get_response(raw_quartiles, par)
  new_iqr <- diff(new_quartiles)
  new_iqr
}

# get_monotonic_penalty <- function(new_scores) {
#   if (any(diff(new_scores$new_score) < 0)) 1000 else 0
# }

get_new_scores <- function(df, par) {
  df %>% 
    mutate(
      new_score = get_response(raw_score, par),
      new_iqr = map2_dbl(raw_score, raw_sem, rescale_error, par)
    )
}

get_iqr_variance <- function(new_scores) {
  var(new_scores$new_iqr)
}

plot_new_iqr_by_score <- function(new_scores) {
  ggplot(new_scores, aes(new_score, new_iqr)) + 
    geom_line() + 
    scale_x_continuous("New score") + 
    scale_y_continuous("New IQR", limits = c(0, NA)) +
    theme(aspect.ratio = 1)
}

eval_par <- function(par, df, plot = FALSE) {
  par <- norm_par(par)
  roughness <- get_roughness(par)
  p1 <- plot_par(par) + ggtitle(sprintf("Roughness = %.2f", roughness))
  new_scores <- get_new_scores(df, par)
  iqr_variance <- get_iqr_variance(new_scores)
  p2 <- plot_new_iqr_by_score(new_scores) + ggtitle(sprintf("IQR SD = %.2f", sqrt(iqr_variance))) 
  if (plot) print(cowplot::plot_grid(p1, p2, nrow = 1))
  cost <- LAMBDA_SQR * roughness + iqr_variance
  message("cost = ", cost)
  cost
}


########################################################

setup_workspace()

RANGE <- c(-10, 10)
INT <- 0.5
BREAKS <- seq(from = RANGE[1], to = RANGE[2], by = INT)
NORDER <- 4
LAMBDA_SQR <- 1e-3

NBASIS <- length(BREAKS) + NORDER - 2L

df_all <- master %>%
  select(raw_score = PIAT.score,
         raw_sem = PIAT.error) %>%
  mutate_all(as.numeric) %>% 
  na.omit()

get_final_sem <- function(raw_score, raw_sem, final_model) {
  N <- 1e4
  raw_samples <- rnorm(N, mean = raw_score, sd = raw_sem)
  transformed_samples <- predict(final_model, newdata = tibble(raw_score = raw_samples))
  sd(transformed_samples)
}

#######

# Original relationship between score and sem
plot(df_all$raw_score, df_all$raw_sem,
     xlab = "Raw score",
     ylab = "Raw SEM")

# Summarised relationship beween score and sem
mod <- gam(raw_sem ~ s(raw_score), data = df_all)
plot(mod)

df <- tibble(
  raw_score = seq(from = -4, to = 4, by = 0.25),
  raw_sem = predict(mod, newdata = tibble(raw_score))
)
plot(df, type = "l", 
     xlab = "Raw score",
     ylab = "Raw SEM")

initial_par <- seq(from = 1, to = 2, length.out = NBASIS)
res <- optim(initial_par, eval_par, df = df, method = "Nelder-Mead", control = list(maxit = 5e3))
par <- norm_par(res$par)
eval_par(par, df, plot = T)

new_scores <- get_new_scores(df, par)

mod_2 <- scam(new_score ~ s(raw_score, bs = "mpi"), data = new_scores)
plot(mod_2)

shift <- as.numeric(predict(mod_2, newdata = tibble(raw_score = 0)))

final_scores <- new_scores %>% 
  mutate(final_score = as.numeric(predict(mod_2, newdata = new_scores %>% select(raw_score))) - shift,
         final_sem = map2_dbl(raw_score, raw_sem, get_final_sem, mod_2))

tibble(
  raw_score = seq(from = -10, to = 10, length.out = 1e2),
  final_score = predict(mod_2, newdata = tibble(raw_score = raw_score)) - shift
) %>% 
  ggplot(aes(raw_score, final_score)) +
  geom_line()

p <- ggplot(final_scores, aes(raw_score, final_score)) +
  geom_line() + 
  scale_x_continuous("Raw score") +
  scale_y_continuous("Final score")
print(p)

q <- ggplot(final_scores, aes(raw_score, new_iqr)) + 
  geom_line() + 
  scale_x_continuous("Raw score") + 
  scale_y_continuous("Final SEM", limits = c(0, NA))
print(q)

# summary <- list(
#   stage_1_mod = norm_par(res),
#   stage_2_mod = mod_2,
#   scores = final_scores
# )