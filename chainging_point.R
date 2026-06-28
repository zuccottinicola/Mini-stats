rm(list = ls())
library(coda)


data <- read.csv("assignment2_changepoint_counts.csv", header = TRUE)
x <- data[,2]
N <- length(x)

# Media mobile
window <- 10
ma <- filter(x, rep(1/window, window), sides = 2)

# Serie temporale
plot(x, type = "l",col  = "steelblue", main = "Time series",  xlab = "i", ylab = "x")


abline(h   = mean(x), col = "red",lwd = 2,lty = 2)

lines(ma,col = "darkred",lwd = 2)

legend("topleft",
       legend = c("Observed", "Global mean", paste0("Moving avg (k=", window, ")")),
       col = c(adjustcolor("steelblue", alpha.f = 0.5), "red", "darkred"),
       lwd = c(1, 2, 2),
       lty = c(1, 2, 1),
       bty = "n",
       cex = 0.85)


hist(x, breaks = 10, main ="Observed counts", col = "lightblue")

alpha <- 1
beta  <- 0.1

nmin <- 5
Mset <- nmin:(N - nmin)

S <- cumsum(x)
SN <- sum(x)

gibbs_cp <- function(x,n_iter = 20000,m_init = 20,lambda1_init = 1,lambda2_init = 5,alpha = 1,beta = 0.1,nmin = 5){
  
  N <- length(x)
  
  Mset <- nmin:(N - nmin)
  
  S <- cumsum(x)
  SN <- sum(x)
  
  
  lambda1 <- numeric(n_iter)
  lambda2 <- numeric(n_iter)
  m_chain <- numeric(n_iter)
  
  
  lambda1[1] <- lambda1_init
  lambda2[1] <- lambda2_init
  m_chain[1] <- m_init

  
  for(t in 2:n_iter){
    
    m_curr <- m_chain[t-1]
    
    shape1 <- alpha + S[m_curr]
    rate1  <- beta + m_curr
    
    lambda1[t] <- rgamma(1,shape = shape1,rate  = rate1)
    
    shape2 <- alpha + (SN - S[m_curr])
    rate2  <- beta + (N - m_curr)
    
    lambda2[t] <- rgamma(1,shape = shape2, rate  = rate2)
    
    log_probs <- numeric(length(Mset))
    
    for(j in seq_along(Mset)){
      
      m <- Mset[j]
      
      S1 <- S[m]
      S2 <- SN - S[m]
      
      log_probs[j] <- S1*log(lambda1[t])-m*lambda1[t]+S2*log(lambda2[t])-(N-m)*lambda2[t]
    }
    
    # Numerical stabilization
    log_probs <- log_probs - max(log_probs)
    
    probs <- exp(log_probs)
    
    probs <- probs / sum(probs)
    
    m_chain[t] <- sample(Mset,size = 1, prob = probs)
  }
  
  return(list(lambda1 = lambda1,lambda2 = lambda2,m= m_chain))
}

set.seed(123)
# 3 chains
chain1 <- gibbs_cp(x,n_iter = 25000,m_init = 10,lambda1_init = 1,lambda2_init = 5)
chain2 <- gibbs_cp(x,n_iter = 25000,m_init = floor(N/2),lambda1_init = 8,lambda2_init = 2)
chain3 <- gibbs_cp(x, n_iter = 25000, m_init = N - 10,lambda1_init = 3, lambda2_init = 9)

burnin <- 5000

post1 <- list(lambda1 = chain1$lambda1[-(1:burnin)],lambda2 = chain1$lambda2[-(1:burnin)],m = chain1$m[-(1:burnin)])
post2 <- list(lambda1 = chain2$lambda1[-(1:burnin)],lambda2 = chain2$lambda2[-(1:burnin)],m = chain2$m[-(1:burnin)])
post3 <- list(lambda1 = chain3$lambda1[-(1:burnin)],lambda2 = chain3$lambda2[-(1:burnin)],m = chain3$m[-(1:burnin)])


par(mfrow = c(3, 3),mar = c(3, 3, 2.5, 1), mgp= c(1.8, 0.5, 0), bg= "white",family = "serif")

make_trace <- function(samples, title, color, burnin = 5000) {
  plot(samples,type = "l",col  = adjustcolor(color, alpha.f = 0.8),
       lwd = 0.8,
       main = title,
       xlab = "Iteration",
       ylab = "",
       cex.main = 1.1,
       cex.axis = 0.85)
  abline(v = burnin, col ="black", lty = 2, lwd = 1.5)}

# Chain 1
make_trace(chain1$lambda1, expression("Chain 1: " * lambda[1]), "steelblue")
make_trace(chain1$lambda2, expression("Chain 1: " * lambda[2]), "firebrick")
make_trace(chain1$m, "Chain 1: m", "darkgreen")

# Chain 2
make_trace(chain2$lambda1, expression("Chain 2: " * lambda[1]), "steelblue")
make_trace(chain2$lambda2, expression("Chain 2: " * lambda[2]), "firebrick")
make_trace(chain2$m, "Chain 2: m", "darkgreen")

# Chain 3
make_trace(chain3$lambda1, expression("Chain 3: " * lambda[1]), "steelblue")
make_trace(chain3$lambda2, expression("Chain 3: " * lambda[2]), "firebrick")
make_trace(chain3$m, "Chain 3: m", "darkgreen")

#unisce chain
lambda1_all <- c(post1$lambda1, post2$lambda1, post3$lambda1)
lambda2_all <- c(post1$lambda2, post2$lambda2, post3$lambda2)
m_all <- c(post1$m, post2$m, post3$m)


par(mfrow = c(1, 3))

hist(lambda1_all, breaks=40, probability=TRUE, col="lightblue", main=expression("Posterior of "*lambda[1]), xlab=expression(lambda[1]))

hist(lambda2_all, breaks=40, probability=TRUE, col="salmon", main=expression("Posterior of "*lambda[2]), xlab=expression(lambda[2]))

hist(m_all, breaks=40, probability=TRUE, col="lightgreen", main="Posterior of m", xlab="m")


par(mfrow = c(1, 3), mar= c(4, 4, 3, 1), bg= "white",family = "serif")

my_acf <- function(x, title, color = "steelblue") {
  acf_obj <- acf(x, plot = FALSE, lag.max = 40)
  ci <- qnorm(0.975) / sqrt(length(x))
  
  plot(acf_obj$lag, acf_obj$acf,
       type     = "h",
       lwd      = 2,
       col      = color,
       ylim     = c(-0.1, 1),
       xlab     = "Lag",
       ylab     = "ACF",
       main     = title,
       cex.main = 1.1,
       cex.axis = 0.85)
  
  abline(h =  ci, lty = 2, col = "red", lwd = 1.5)
  abline(h = -ci, lty = 2, col = "red", lwd = 1.5)
  abline(h =  0,  lty = 1, col = "black", lwd = 0.8)
  
  points(acf_obj$lag, acf_obj$acf, pch = 19, cex = 0.6, col = color)
}

my_acf(lambda1_all, expression("ACF of " * lambda[1]), color = "steelblue")
my_acf(lambda2_all, expression("ACF of " * lambda[2]), color = "firebrick")
my_acf(m_all, "ACF of m", color = "darkgreen")

ess_lambda1 <- effectiveSize(lambda1_all)
ess_lambda2 <- effectiveSize(lambda2_all)
ess_m <- effectiveSize(m_all)

cat("\nEffective Sample Sizes:\n")

cat("ESS lambda1 =", ess_lambda1, "\n")
cat("ESS lambda2 =", ess_lambda2, "\n")
cat("ESS m =", ess_m, "\n")


posterior_summary <- function(samples){
  
  density_est <- density(samples)
  MAP <- density_est$x[which.max(density_est$y)]
  
  c(Mean= mean(samples), Median = median(samples), SD= sd(samples), MAP= MAP, CI_low = quantile(samples, 0.025), CI_up = quantile(samples, 0.975))
}

summary_lambda1 <- posterior_summary(lambda1_all)
summary_lambda2 <- posterior_summary(lambda2_all)
summary_m <- posterior_summary(m_all)

cat("\nPosterior summary lambda1:\n")
print(summary_lambda1)

cat("\nPosterior summary lambda2:\n")
print(summary_lambda2)

cat("\nPosterior summary m:\n")
print(summary_m)


prob <- mean(lambda2_all > lambda1_all)
cat("\nP(lambda2 > lambda1 | D) =", prob, "\n")



x_new <- rpois(length(lambda2_all),lambda = lambda2_all)
x_bar <- mean(x)

prob_large <- mean(x_new > 2 * x_bar)

cat("\nPosterior predictive mean = ",mean(x_new),"\n")

cat("\nP(x_new > 2*x_bar | D) = ",  prob_large,"\n")
par(mfrow = c(1, 1))

hist(x_new,breaks = 30, probability = TRUE,col = "lightyellow",main = "Posterior predictive distribution",xlab = expression(x[new]))

abline(v = x_bar,col = "darkgreen", lty = 2, lwd = 1.5)

abline(v = 2*x_bar,  col = "red",lty = 2, lwd = 1.5)

legend("topright",legend = c(expression(bar(x)), expression(2*bar(x))), col = c("darkgreen","red"),lwd = 2)

log_marginal_M0 <- function(x, alpha, beta){
  
  N  <- length(x)
  SN <- sum(x)
  
  lgamma(alpha + SN) -(alpha + SN)*log(beta + N) +alpha*log(beta)-lgamma(alpha)
}

log_marginal_M1 <- function(x, alpha, beta, nmin){
  
  N  <- length(x)
  S  <- cumsum(x)
  Mset <- nmin:(N - nmin)
  log_terms <- numeric(length(Mset))
  
  for(j in seq_along(Mset)){
    
    m <- Mset[j]
    
    S1 <- S[m]
    S2 <- S[N] - S[m]
    
    log_terms[j] <- lgamma(alpha+S1)-(alpha+S1)*log(beta+m)+lgamma(alpha+S2)-(alpha+S2)*log(beta+N-m)
  }
  
  max_log <- max(log_terms)
  
  max_log + log(sum(exp(log_terms - max_log)))
}

logBF10 <- log_marginal_M1(x, alpha, beta, nmin)-log_marginal_M0(x, alpha, beta)

BF10 <- exp(logBF10)

cat("\nBayes Factor = ", BF10, "\n")
