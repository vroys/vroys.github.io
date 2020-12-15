#######################################################
##R codes for the short course on Monte Carlo methods
##and Bayesian statistics given at Chennai Mathematical
##Institute, June 10, 2019
##By Vivekananda Roy (vroy@iastate.edu)
#######################################################
## Monte Carlo Toy Examples
#######################################################

set.seed(3)

n <- 1000
x <- rnorm(n, me=1)
y <- sqrt(2*pi)*x
est <- mean(y)
est
mcse <- sd(y) / sqrt(n)
interval <- est + c(-1,1)*1.96*mcse
interval


y <- sqrt(2*pi)*x*sin(x)
est <- mean(y)
est
mcse <- sd(y) / sqrt(n)
interval <- est + c(-1,1)*1.96*mcse
interval
#sqrt(2*pi)*exp(-1/2)*(sin(1)+cos(1))


n <- 1000
x <- rexp(n, rate=.5)
y <- x^2
est <- mean(y)
est
mcse <- sd(y) / sqrt(n)
interval <- est + c(-1,1)*1.96*mcse
interval

y <- x^2/log(x+2)
est <- mean(y)
est
mcse <- sd(y) / sqrt(n)
interval <- est + c(-1,1)*1.96*mcse
interval

######################################################
## MH RW sampler for normal
######################################################

set.seed(3)
library(mcmcse)

n_iterations = 10000
sigma=2.4
log_pi = function(y) {
dnorm(y,log=TRUE)
}
current = -1 # Initial value
samps = rep(NA,n_iterations)
for (i in 1:n_iterations) {
proposed = rnorm(1, current, sigma)
logr = log_pi(proposed)-log_pi(current)
if (log(runif(1)) < logr) current = proposed
samps[i] = current
}

length(unique(samps))/n_iterations # acceptance rate

ts.plot(samps[1:1000])


mcse(samps)

plot(acf(samps))

######################################################
## MH independence sampler for normal
######################################################

set.seed(3)
library(mcmcse)

n_iterations = 10000
sigma=2.4
log_pi = function(y) {
dnorm(y,log=TRUE)
}
log_f = function(y) {
dnorm((y/sigma),log=TRUE)-log(sigma)
}
current = 0.5 # Initial value
samps = rep(NA,n_iterations)
for (i in 1:n_iterations) {
proposed = rnorm(1,me=0,sd=sigma)
logr = log_pi(proposed)+log_f(current)-log_pi(current)-log_f(proposed)
if (log(runif(1)) < logr) current = proposed
samps[i] = current
}

length(unique(samps))/n_iterations # acceptance rate

ts.plot(samps[1:1000])

mcse(samps)

plot(acf(samps))

######################################################
## Gibbs sampler for bivariate normal
######################################################
set.seed(3)
library(mcmcse)
gibbs_bivariate_normal = function(samps_start, n_iterations, rho) {
    samps = matrix(samps_start, nrow=n_iterations, ncol=2, byrow=TRUE)
v = sqrt(1-rho^2)
for (i in 2:n_iterations) {
samps[i,1] = rnorm(1, rho*samps[i-1,2], v)
samps[i,2] = rnorm(1, rho*samps[i ,1], v)
}
return(samps)
}
samps = gibbs_bivariate_normal(samps_start<-c(-3,3), n_iterations<-1000, rho<-0.9)

ts.plot(samps[,1])
ts.plot(samps[,2])
plot(acf(samps[,2]))
apply(samps, 2, mcse)
apply(samps, 2, quantile, probs=c(0.025,0.975))
cor(samps[,1],samps[,2])
##########################################################
## Gibbs sampler for Normal Regression with Jeffreys prior
##########################################################
set.seed(3)

library(MCMCpack)
library(mcmcse)
data<-read.csv("stock_treasury.csv")
# Risk Free Rate is in percentage and annualised.
# So the following conversion is required.
Rf<-data$UST_Yr_1/(100*250)
plot(ts(Rf),ylab="US Treasury 1 Year Yield")
n<-nrow(data)
## Compute log-return
ln_rt_snp500<-diff(log(data$SnP500))-Rf[2:n]
ln_rt_ibm<-diff(log(data$IBM_AdjClose))-Rf[2:n]
ln_rt_apple<-diff(log(data$Apple_AdjClose))-Rf[2:n]
ln_rt_msft<-diff(log(data$MSFT_AdjClose))-Rf[2:n]
ln_rt_intel<-diff(log(data$Intel_AdjClose))-Rf[2:n]

y = ln_rt_ibm
n_obs =length(y)
X = cbind(rep(1, n_obs), ln_rt_snp500) #include an intercept
XtX = t(X) %*% X
n_params = 2
n_obsby2=n_obs/2

beta_hat = solve(XtX, t(X) %*% y) # compute this beforehand
XtXi = solve(XtX)

beta = c(0,0) # starting value
                                        #beta =beta_hat

n_iterations = 5000 #number of MCMC iterations


beta_out = matrix(data=NA, nrow=n_iterations, ncol=n_params)
sigma_out = matrix(data = NA, nrow = n_iterations, ncol=1)

for (i in 1:n_iterations){
      ymxbeta= (y - X %*% beta)
     sigma2 = rinvgamma(1, n_obsby2, t(ymxbeta) %*% ymxbeta * .5 ) # draw from sigma2 given beta
      beta = mvrnorm(n=1, beta_hat, sigma2 * XtXi) # draw from beta given sigma2

      ##save iterations
  sigma_out[i,] = sigma2
  beta_out[i,] = beta
}

ts.plot(sigma_out)
ts.plot(beta_out[,2])
apply(beta_out,2,mean)
mean(sigma_out)
apply(beta_out, 2, mcse)
mcse(sigma_out)
apply(beta_out, 2, quantile, probs=c(0.025,0.975))
quantile(sigma_out, probs=c(0.025,0.975))

