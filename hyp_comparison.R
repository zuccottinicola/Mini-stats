data <- scan("dataset.txt")
N <- length(data)
x_bar <- mean(data)
sum_x <- sum(data)

summary(data)
hist(data, breaks = 20, main = "Observed counts", col = "lightblue")


# Prior hyperparameters
alpha <- 1
beta  <- 1

# Posterior parameters
alpha_post <- alpha + sum_x
beta_post  <- beta + N
lambda_map <- (alpha_post - 1) /(beta_post)
lambda_map

# Posterior mean and variance
lambda_mean <- alpha_post / beta_post
lambda_var  <- alpha_post / (beta_post^2)
lambda_sd   <- sqrt(lambda_var)

lambda_mean
lambda_sd

# Plot posterior
curve(dgamma(x, shape = alpha_post, rate = beta_post),
      from = 0, to = lambda_mean + 4*lambda_sd,
      main = "Posterior of lambda", ylab = "Density")

abline(v = lambda_map, col = "red", lwd = 2)


#ex 2

threshold <- 2 * x_bar

k <- floor(threshold)

prob_exceed <- 1 - pnbinom(k,size = alpha_post, prob = beta_post / (beta_post + 1))

prob_exceed


# ex 3

log_likelihood_nb <- function(params, data) {
  r <- params[1]
  pi <- params[2]
  
  if (r <= 0 || pi <= 0 || pi >= 1) return(-Inf)
  
  sum(dnbinom(data, size = r, prob = pi, log = TRUE))
}

log_prior_nb <- function(params) {
  r <- params[1]
  pi <- params[2]
  
  if (r <= 0 || pi <= 0 || pi >= 1) return(-Inf)
  
  a <- 1; b <- 1
  c <- 1; d <- 1
  
  dgamma(r, shape = a, rate = b, log = TRUE) + dbeta(pi, shape1 = c, shape2 = d, log = TRUE)
}

log_posterior_nb <- function(params, data) {
  log_likelihood_nb(params, data) + log_prior_nb(params)
}

init <- c(r = 1, pi = 0.5)

fit <- optim(init,fn = function(p) -log_posterior_nb(p, data),hessian = TRUE)

r_map  <- fit$par[1]
pi_map <- fit$par[2]

r_map
pi_map
cov_matrix <- solve(fit$hessian)
se <- sqrt(diag(cov_matrix))

se_r  <- se[1]
se_pi <- se[2]

se_r
se_pi
k_param <- 2

log_evidence_nb <- -fit$value +(k_param/2)*log(2*pi) + 0.5 * log(det(cov_matrix))

log_evidence_nb
log_evidence_pois <- alpha * log(beta) -lgamma(alpha) + lgamma(alpha_post) -(alpha_post) * log(beta_post)-sum(lgamma(data + 1))

log_evidence_pois
log_BF <- log_evidence_nb - log_evidence_pois
BF <- exp(log_BF)

BF



#ex 4
if (BF > 1) {
  cat("Negative Binomial favored\n")
} else {
  cat("Poisson favored\n")
}


x_vals <- 0:max(data)

nb_pmf <- dnbinom(x_vals, size = r_map, prob = pi_map)
hist(data,
     probability = TRUE,
     breaks = 20,
     col  = "grey",
     main = "Fit: Poisson vs Negative Binomial",
     xlim = c(0, max(x_vals)),                          
     ylim = c(0, max(c(nb_pmf, pois_pmf)) * 1.1)) 
#hist(data, probability = TRUE, breaks = 20, col = "grey", main = "Fit: Poisson vs Negative Binomial")

points(x_vals, nb_pmf, col = "red", pch = 16)
lines(x_vals, nb_pmf, col = "red", lwd = 2)
pois_pmf <- dpois(x_vals, lambda_map)

points(x_vals, pois_pmf, col = "blue", pch = 16)
lines(x_vals, pois_pmf, col = "blue", lwd = 2)

legend("topright",legend = c("Poisson", "NegBin"),col = c("blue", "red"),lwd = 2)

