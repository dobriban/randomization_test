# Experiment with mean and permutation version
# Generating data

#options: "t" or "mean"
#stat_type = "mean"
stat_type = "t"

set.seed(2)

alpha <- 0.05
p <- 10


num_mu <- 20
mus <- seq(0, 3, 3 / (num_mu - 1))
# num_mu <- 1
# mus <- 0


pow_t = rep(0, num_mu)
pow_perm =  rep(0, num_mu)

num_sim <- 1000#1000
num_perm <- 1000

#Permutation test fn
one.test <- function(ind, y) {
  xstar <- sample(ind)
  res <- stat(y[xstar == 1], y[xstar == 0])
  return(res)
}

if (stat_type == "t") {
  stat <- function(x1, x2) {
    p <- length(x1)
    tt <- t.test(x1, x2, var.equal = T)
    s_det <- tt$statistic
    p_det <- tt$p.value
    return(c(s_det, p_det))
  }
}
if (stat_type == "mean") {
  stat <- function(x1, x2) {
    p <- length(X1)
    s_det <- (mean(x1) - mean(x2)) / sqrt(2 / p)
    p_det <- 1 - pnorm(s_det)
    return(c(s_det, p_det))
  }
}


for (i in 1:num_mu) {
  mu <- mus[i]
  
  rej_t <- rep(0, num_sim)
  rej_perm <- rep(0, num_sim)
  
  for (j in 1:num_sim) {
    X1 <- rnorm(p) + mu * rep(1, p)
    X2 <- rnorm(p)
    s <- stat(X1, X2)
    original_stat <- s[1]
    rej_t[j] <- 1 * (s[2] < alpha)
    
    #T-test
    #this code partly inspired by https://faculty.washington.edu/kenrice/sisg/SISG-08-06.pdf
    X <- c(X1, X2)
    ind <- c(rep(1, p), rep(0, p))
    perm_stats <- replicate(num_perm, one.test(ind, X))
    pval <- (sum(perm_stats >= original_stat) + 1) / (num_perm + 1)
    
    
    rej_perm[j] <- 1 * (pval < alpha)
  }
  
  pow_t[i] = mean(rej_t)
  pow_perm[i] = mean(rej_perm)
  
}


draw <- 0
if (draw == 1) {
  hist(pvals, 100)
}

########
if (stat_type == "t") {
  png(file = "t_exp.png",
      width = 600,
      height = 350)
}

if (stat_type == "mean") {
  png(file = "mean_exp.png",
      width = 600,
      height = 350)
}


plot(
  x = mus,
  y = pow_t,
  type = "l",
  lty = 1,
  ylim = c(0, 1),
  axes = T,
  bty = "n",
  xaxs = "i",
  yaxs = "i",
  xlab = "mu",
  ylab = "Power"
)

# plot dashed line
lines(x = mus, y = pow_perm, lty = 2)

# add legend
par(xpd = TRUE)
legend(
  x = 1.5,
  y = 0.2,
  legend = c("Deterministic", "Permutation"),
  lty = 1:2,
  box.lty = 0,
  ncol = 2
)
dev.off()
