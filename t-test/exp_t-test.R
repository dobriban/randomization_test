# Generating data
set.seed(2)

alpha <- 0.05
p <- 10
num_mu <- 20
mus <- seq(0, 3, 3 / (num_mu-1))

pow_t = rep(0, num_mu)
pow_perm =  rep(0, num_mu)

num_sim <- 100#1000
num_perm <- 1000

#Permutation test fn
permutation.test <- function(X1, X2, num_perm) {
  permuted_stat = c()
  result = 0
  for (i in 1:num_perm) {
    X <- c(X1, X2)
    Xn <- split(sample(X), rep(1:2, c(p, p)))
    tt <- t.test(Xn[[1]], Xn[[2]], var.equal = T)
    permuted_stat[i] = tt$statistic
  }
  pval = sum(permuted_stat > original_stat) / (num_perm + 1)
  return(list(pval, permuted_stat))
}

one.test <- function(ind,y) {
  xstar<-sample(ind)
  t.test(y[xstar==1], y[xstar==0], var.equal = T)$statistic
}

for (i in 1:num_mu) {
  mu <- mus[i]
  
  rej_t <- rep(0, num_sim)
  rej_perm <- rep(0, num_sim)

  for (j in 1:num_sim) {
    X1 <- rnorm(p) + mu * rep(1, p)
    X2 <- rnorm(p)
    tt <- t.test(X1, X2, var.equal = T)
    original_stat <- tt$statistic
    rej_t[j] <- 1 * (tt$p.value < alpha)    
    
    #T-test: version 1
    # perm_tt <- permutation.test(X1, X2, num_perm)
    # pval <- perm_tt[[1]]
    
    #T-test: version 2
    #this code partly inspired by https://faculty.washington.edu/kenrice/sisg/SISG-08-06.pdf
    X <- c(X1,X2)
    ind <- c(rep(1,p),rep(0,p))
    perm_stats <- replicate(num_perm, one.test(ind, X))
    pval <- mean(perm_stats > original_stat)
    
    
    rej_perm[j] <- 1 * (pval < alpha)
  }
  
    pow_t[i] = mean(rej_t)
    pow_perm[i] = mean(rej_perm)
  
}

########
#png(file="t_exp.png", width=600, height=350)
png(file="t_exp_2.png", width=600, height=350)
plot(x=mus, y=pow_t, type="l", lty=1, ylim=c(0,1),
     axes=T, bty="n", xaxs="i", yaxs="i",
     xlab="mu", ylab="Power")

# plot dashed line
lines(x=mus, y=pow_perm, lty=2)

# add axes
#axis(side=1, labels=, at=1:num_mu)
#axis(side=2, at=seq(0,1,10), las=1)

# add legend
par(xpd=TRUE)
legend(x=1.5, y=0.2, legend=c("Deterministic", "Permutation"), lty=1:2, box.lty=0, ncol=2)
dev.off()





