library(rjags)
library(coda)
library(ggplot2)
library(patchwork)
library(gridExtra)

theme_set(
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "gray40", size = 11),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
)

COL_B0  <- "#2166AC"   
COL_B1  <- "#D6604D"   
COL_FIT <- "#1B7837"   
COL_PI  <- "#762A83"  
COL_DAT <- "gray35"    

set.seed(42)

data_df <- read.csv("logistic_data.csv")
y <- data_df$y
x <- data_df$x
N <- length(y)

x_mean <- mean(x)
x_sd   <- sd(x)
x_std  <- as.numeric(scale(x))

cat("Dataset: N =", N,
    "| successi:", sum(y),
    "| media(y) =", round(mean(y), 3),
    "| x in [", round(min(x), 2), ",", round(max(x), 2), "]\n")


p1b <- ggplot(data_df, aes(x = x, fill = factor(y), colour = factor(y))) +
  geom_density(alpha = 0.35, linewidth = 0.9) +
  scale_fill_manual(values   = c("0" = COL_DAT, "1" = COL_B0),
                    labels   = c("y = 0", "y = 1"), name = NULL) +
  scale_colour_manual(values = c("0" = COL_DAT, "1" = COL_B0),
                      labels = c("y = 0", "y = 1"), name = NULL) +
  labs(title = "Distribution per class",
       x = "x", y = "Density")
p1b

#1 and 2

model_M1 <- "
model {
  for (i in 1:N) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <- beta0 + beta1 * x[i]
  }
  beta0 ~ dnorm(0, 1/9)
  beta1 ~ dnorm(0, 1/9)
}
"

model_M0 <- "
model {
  for (i in 1:N) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <- beta0
  }
  beta0 ~ dnorm(0, 1/9)
}
"

data_M1 <- list(N = N, x = x_std, y = y)
data_M0 <- list(N = N, y = y)


inits_M1 <- list(
  list(beta0 = -5, beta1 = -5),
  list(beta0 = -5, beta1 =  5),
  list(beta0 =  5, beta1 = -5),
  list(beta0 =  5, beta1 =  5)
)
inits_M0 <- list(
  list(beta0 = -5), list(beta0 = -1),
  list(beta0 =  1), list(beta0 =  5)
)


jmod_M1 <- jags.model(textConnection(model_M1),data = data_M1, inits = inits_M1,n.chains = 4, quiet = FALSE)
update(jmod_M1, 1000) #burn in
samp_M1 <- coda.samples(jmod_M1, variable.names = c("beta0", "beta1"), n.iter = 20000)

jmod_M0 <- jags.model(textConnection(model_M0),data = data_M0, inits = inits_M0,n.chains = 4, quiet = TRUE)
update(jmod_M0, 1000)
samp_M0 <- coda.samples(jmod_M0, variable.names = "beta0", n.iter = 20000)



gr  <- gelman.diag(samp_M1)
ess <- effectiveSize(samp_M1)

cat("\nGelman-Rubin R-hat:\n");   print(gr)
cat("\nEffective Sample Size:\n"); print(ess)
cat("\nAutocorrelazione (lag 1,5,10,50):\n")
print(autocorr.diag(samp_M1, lags = c(1, 5, 10, 50)))
cat("\nSintesi posteriore:\n");   print(summary(samp_M1))

draws_list <- lapply(seq_along(samp_M1), function(k) {
  df <- as.data.frame(samp_M1[[k]])
  df$chain <- k
  df$iter  <- seq_len(nrow(df))
  df
})
draws_long <- do.call(rbind, draws_list)

p_tr_b0 <- ggplot(draws_long, aes(iter, beta0, colour = factor(chain))) +
  geom_line(alpha = .7, linewidth = .25) +
  geom_vline(xintercept = 1000, linetype = 2, linewidth = .8,colour = "black") +
  labs(title = expression("Trace plot of " * beta[0]), x = "Iteration", y = expression(beta[0]))

p_tr_b1 <- ggplot(draws_long, aes(iter, beta1, colour = factor(chain))) +
  geom_line(alpha = .7, linewidth = .25) +
  geom_vline(xintercept = 1000, linetype = 2, linewidth = .8,colour = "black") +
  labs(title = expression("Trace plot of " * beta[1]), x = "Iteration", y = expression(beta[1]))

(p_tr_b0|p_tr_b1)

make_acf_plot <- function(chain_vec, col, pname) {
  a   <- acf(chain_vec, lag.max = 60, plot = FALSE)
  ci  <- qnorm(0.975) / sqrt(length(chain_vec))
  adf <- data.frame(lag = a$lag[-1], acf = a$acf[-1])
  ggplot(adf, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, colour = "gray50") +
    geom_hline(yintercept = c(-ci, ci),linetype = "dashed", colour = "steelblue") +
    geom_segment(aes(xend = lag, yend = 0), colour = col, linewidth = 0.8) +
    ylim(-0.2, 1) +
    labs(title = bquote("ACF of " ~ .(pname)),x = "lag", y = "ACF")
}
draws  <- as.matrix(samp_M1)

beta0_all <- as.numeric(draws[, "beta0"])
beta1_all <- as.numeric(draws[, "beta1"])
p_acf_b0 <- make_acf_plot(beta0_all,COL_B0, quote(beta[0]))
p_acf_b1 <- make_acf_plot(beta1_all,COL_B1, quote(beta[1]))

(p_acf_b0 | p_acf_b1)

# 3
beta0  <- as.numeric(draws[, "beta0"])
beta1  <- as.numeric(draws[, "beta1"])

post_tab <- data.frame(
  Parametro = c("beta0", "beta1"),
  Media     = round(c(mean(beta0), mean(beta1)), 4),
  Mediana   = round(c(median(beta0), median(beta1)),  4),
  SD        = round(c(sd(beta0), sd(beta1)), 4),
  CI_2.5    = round(c(quantile(beta0, 0.025), quantile(beta1, 0.025)), 4),
  CI_97.5   = round(c(quantile(beta0, 0.975), quantile(beta1, 0.975)), 4)
)

print(post_tab, row.names = FALSE)
cat("\nCorrelazione posteriore (beta0, beta1):", round(cor(beta0, beta1), 4), "\n")


post_df <- data.frame(beta0 = beta0, beta1 = beta1)
make_marginal <- function(samples, col, pname, ci_lo, ci_hi, mu) {
  
  df <- data.frame(x = samples)
  p <- ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "lightblue", colour = "white") +
    geom_density(linewidth = 1.2) +
    geom_vline(xintercept = mu, colour = col, linewidth = 1.1) +
    geom_vline(xintercept = c(ci_lo, ci_hi), colour = col, linetype = "dashed", linewidth = 0.8) +
    labs(title = bquote("Posterior of " ~ .(pname)),x = pname, y = "Density")
  if (grepl("beta1", deparse(pname))) {
    p <- p + geom_vline(xintercept = 0,colour = "black",linetype = "dotted",linewidth = 0.9)
  }
  
  p
}


p_m0 <- make_marginal(beta0, COL_B0, quote(beta[0]),quantile(beta0, .025), quantile(beta0, .975), mean(beta0))
p_m1 <- make_marginal(beta1, COL_B1, quote(beta[1]), quantile(beta1, .025), quantile(beta1, .975), mean(beta1))
(p_m0|p_m1)


p_joint <- ggplot(post_df, aes(x = beta0, y = beta1)) +
  geom_point(alpha = 0.05, size = 0.4, colour = "steelblue") +
  geom_density_2d(colour = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = mean(beta0), linetype = "dashed", colour = COL_B0) +
  geom_hline(yintercept = mean(beta1), linetype = "dashed", colour = COL_B1) +
  labs(title = expression("Joint posterior: " * beta[0] * " vs " * beta[1]), x = expression(beta[0]),y = expression(beta[1]))

p_joint
# 4

x_new_std  <- 0
x_new_orig <- x_new_std * x_sd + x_mean

pi_new <- plogis(beta0 + beta1 * x_new_std)
y_new  <- rbinom(length(pi_new), size = 1, prob = pi_new)

cat("x_new (std) =", x_new_std, " => x_new (originale) =", round(x_new_orig, 3), "\n")
cat("E[pi_new | D]    =", round(mean(pi_new), 4), "\n")
cat("Mediana pi_new   =", round(median(pi_new), 4), "\n")
cat("CI 95% pi_new    = [", round(quantile(pi_new, .025), 4),",", round(quantile(pi_new, .975), 4), "]\n")
cat("P(y_new=1|x,D)   =", round(mean(y_new), 4), "\n")

pi_stats <- c(media = mean(pi_new),lo95 = as.numeric(quantile(pi_new, .025)),hi95 = as.numeric(quantile(pi_new, .975)), mediana = median(pi_new))
fig4 <- ggplot(data.frame(pi = pi_new), aes(x = pi)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = COL_PI, alpha = 0.3, colour = NA) +
  geom_density(colour = COL_PI, linewidth = 1.0) +
  geom_vline(xintercept = pi_stats["media"], colour = COL_PI, linewidth = 1.0, linetype = "dashed") +
  geom_vline(xintercept = c(pi_stats["lo95"], pi_stats["hi95"]), colour = COL_PI, linewidth = 0.5, linetype = "dashed") +
  labs(title= bquote("Predictive distribution of " ~pi[new]),x = bquote(pi[new]),y = "Density")
fig4
# 5


jmod_M1_dic <- jags.model(textConnection(model_M1),data = data_M1, inits = inits_M1, n.chains = 4, quiet = TRUE)
update(jmod_M1_dic, 1000)
dic_M1 <- dic.samples(jmod_M1_dic, n.iter = 10000, type = "pD")

jmod_M0_dic <- jags.model(textConnection(model_M0),data = data_M0, inits = inits_M0,n.chains = 4, quiet = TRUE)
update(jmod_M0_dic, 1000)
dic_M0 <- dic.samples(jmod_M0_dic, n.iter = 10000, type = "pD")

cat("\nDIC M1:\n"); print(dic_M1)
cat("\nDIC M0:\n"); print(dic_M0)

dd <- diffdic(dic_M1, dic_M0)   
delta_DIC <- as.numeric(dd)

cat("\nDelta DIC (M1 - M0) =", round(delta_DIC,2), "\n")

    