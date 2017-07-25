library(rstanarm)
library(ggplot2)
library(bayesplot)

# Compute Bayesian R-squared for linear and
# generalized linear models.
#
# @param fit A fitted model object returned by stan_glm.
# @return A vector of R-squared values with length equal to
#      the number of posterior draws.
#
bayes_R2 <- function(fit) {
  y <- get_y(fit)
  ypred <- posterior_linpred(fit, transform = TRUE)
  if (family(fit)$family == "binomial" && NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1]
    ypred <- ypred %*% diag(trials)
  }
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}


# Example -----------------------------------------------------------------

x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3
xy <- data.frame(x,y)

## Lsq fit
fit_lm <- lm(y ~ x, data = xy)
ols_coef <- coef(fit_lm)

yhat <- ols_coef[1] + ols_coef[2] * x
r <- y - yhat
rsq_1 <- var(yhat)/(var(y))
rsq_2 <- var(yhat)/(var(yhat) + var(r))

print(c(rsq_1, rsq_2))

## Bayes fit
fit_stan <- stan_glm(
  y ~ x,
  data = xy,
  prior_intercept = normal(0, 0.2, autoscale = FALSE),
  prior = normal(1, 0.2, autoscale = FALSE)
)
posterior <- as.matrix(fit_stan, pars = c("(Intercept)", "x"))
post_means <- colMeans(posterior)

print(median(bayes_R2(fit_stan)))


# Figures -----------------------------------------------------------------
theme_set(bayesplot::theme_default(base_family = "sans"))

fig_1a <-
  ggplot(xy, aes(x, y)) +
  geom_point() +
  geom_abline( # ols regression line
    intercept = ols_coef[1],
    slope = ols_coef[2],
    size = 0.5
  ) +
  geom_abline( # prior regression line
    intercept = 0,
    slope = 1,
    color = "blue",
    linetype = 2,
    size = 0.5
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "blue3",
    size = 0.5
  ) +
  geom_point(
    aes(y = post_means[1] + post_means[2] * x),
    color = "blue"
  ) +
  coord_fixed(xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25)) +
  annotate(
    geom = "text",
    x = c(-1.7, -.85, 1.1),
    y = c(-0.63, -1.85, 1),
    label = c(
      "Least-squares\nfit",
      "Prior regression line",
      "Posterior mean fit"
    ),
    color = c("black", "blue", "blue3"),
    size = 3.8
  ) +
  ggtitle("Least squares and Bayes fits") +
  theme(plot.title = element_text(face = "bold"))

plot(fig_1a)
ggsave("fig/fig_1a.pdf", width = 4.5, height = 4.5)


# take 20 posterior draws
keep <- sample(nrow(posterior), 20)
samp <- posterior[keep, ]

fig_1b <-
  ggplot(xy, aes(x, y)) +
  geom_point() +
  geom_abline( # 20 posterior draws of the regression line
    intercept = samp[, 1],
    slope = samp[, 2],
    color = "blue",
    size = 0.25,
    alpha = 1/3
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "blue3",
    size = 1
  ) +
  coord_fixed(xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25)) +
  ggtitle("Bayes posterior simulations") +
  theme(plot.title = element_text(face = "bold"))

plot(fig_1b)
ggsave("fig/fig_1b.pdf", width = 4.5, height = 4.5)

