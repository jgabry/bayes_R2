# File contents: 
#   * bayes_R2 function
#   * code for fitting the models for the example in the paper
#   * code for producing the plots in the paper (both base graphics and ggplot2 code)


# Compute Bayesian R-squared for stan_glm models
#
# NOTE: If you have a rstanarm v2.15.4 or higher then the bayes_R2
# function is included ---see help("bayes_R2", "rstanarm")--- and
# you do not need to use this simplified version here. The version
# in the rstanarm package is also more compatible with stan_glmer
# models than this version.
#
# @param fit A fitted model object returned by stan_glm.
# @return A vector of R-squared values with length equal to
#      the number of posterior draws.
#
bayes_R2 <- function(fit) {
  y <- rstanarm::get_y(fit)
  ypred <- rstanarm::posterior_linpred(fit, transform = TRUE)
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

library("rstanarm")

x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3
xy <- data.frame(x,y)

## Least squares fit
fit <- lm(y ~ x, data = xy)
ols_coef <- coef(fit)

yhat <- ols_coef[1] + ols_coef[2] * x
r <- y - yhat
rsq_1 <- var(yhat) / (var(y))
rsq_2 <- var(yhat) / (var(yhat) + var(r))

print(c(rsq_1, rsq_2))

## Bayes fit with strong priors
fit_bayes <- stan_glm(
  y ~ x,
  data = xy,
  prior_intercept = normal(0, 0.2, autoscale = FALSE),
  prior = normal(1, 0.2, autoscale = FALSE)
)
posterior <- as.matrix(fit_bayes, pars = c("(Intercept)", "x"))
post_means <- colMeans(posterior)

rsq_bayes <- bayes_R2(fit_bayes)
print(median(rsq_bayes))


# Figures -----------------------------------------------------------------

# The first section of code below creates plots using base R graphics. 
# Below that there is code to produce the plots using ggplot2.

# take a sample of 20 of the posterior draws
keep <- sample(nrow(posterior), 20)
samp_20_draws <- posterior[keep, ]

# base graphics version

pdf("fig/rsquared1a.pdf", height=4, width=5)
par(mar=c(3,3,1,1), mgp=c(1.7,.5,0), tck=-.01)
plot(
  x, y,
  ylim = range(x),
  xlab = "x",
  ylab = "y",
  main = "Least squares and Bayes fits",
  bty = "l",
  pch = 20
)
abline(coef(fit)[1], coef(fit)[2], col = "black")
text(-1.6,-.7, "Least-squares\nfit", cex = .9)
abline(0, 1, col = "blue", lty = 2)
text(-1, -1.8, "(Prior regression line)", col = "blue", cex = .9)
abline(coef(fit_bayes)[1], coef(fit_bayes)[2], col = "blue")
text(1.4, 1.2, "Posterior mean fit", col = "blue", cex = .9)
points(
  x,
  coef(fit_bayes)[1] + coef(fit_bayes)[2] * x,
  pch = 20,
  col = "blue"
)
dev.off()

## I used color-hex.com

pdf("fig/rsquared1b.pdf", height=4, width=5)
par(mar=c(3,3,1,1), mgp=c(1.7,.5,0), tck=-.01)
plot(
  x, y,
  ylim = range(x),
  xlab = "x",
  ylab = "y",
  bty = "l",
  pch = 20,
  main = "Bayes posterior simulations"
)
for (s in 1:nrow(samp_20_draws)) {
  abline(samp_20_draws[s, 1], samp_20_draws[s, 2], col = "#9497eb")
}
abline(
  coef(fit_bayes)[1],
  coef(fit_bayes)[2],
  col = "#1c35c4",
  lwd = 2
)
points(x, y, pch = 20, col = "black")
dev.off()



# ggplot version

library("ggplot2")
library("bayesplot")

theme_set(bayesplot::theme_default(base_family = "sans"))
theme_update(
  plot.title = element_text(face = "bold", hjust = 0.5), 
  axis.text = element_text(size = rel(1.1))
)

fig_1a <-
  ggplot(xy, aes(x, y)) +
  geom_point() +
  geom_abline( # ols regression line
    intercept = ols_coef[1],
    slope = ols_coef[2],
    size = 0.4
  ) +
  geom_abline( # prior regression line
    intercept = 0,
    slope = 1,
    color = "blue",
    linetype = 2,
    size = 0.3
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "blue",
    size = 0.4
  ) +
  geom_point(
    aes(y = post_means[1] + post_means[2] * x),
    color = "blue"
  ) +
  annotate(
    geom = "text",
    x = c(-1.6, -1, 1.4),
    y = c(-0.7, -1.8, 1.2),
    label = c(
      "Least-squares\nfit",
      "(Prior regression line)",
      "Posterior mean fit"
    ),
    color = c("black", "blue", "blue"),
    size = 3.8
  ) +
  ylim(range(x)) + 
  ggtitle("Least squares and Bayes fits")

plot(fig_1a)
ggsave("fig/rsquared1a-gg.pdf", width = 5, height = 4)


fig_1b <-
  ggplot(xy, aes(x, y)) +
  geom_abline( # 20 posterior draws of the regression line
    intercept = samp_20_draws[, 1],
    slope = samp_20_draws[, 2],
    color = "#9497eb",
    size = 0.25
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "#1c35c4",
    size = 1
  ) +
  geom_point() +
  ylim(range(x)) + 
  ggtitle("Bayes posterior simulations")

plot(fig_1b)
ggsave("fig/rsquared1b-gg.pdf", width = 5, height = 4)

