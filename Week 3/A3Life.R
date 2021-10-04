setwd(paste0(getwd(),"/Week 3/"))
rm(list = ls())

# Question 4
gen.Sample <- function(n, a, b, c) {
  if (any(a<0, b<0, c<1)) stop("Invalid parameters")
  lifetimes <- log(1+rexp(n)*log(c)/b)/log(c) ## Gompertz sample; see Q3
  if (a>0) lifetimes <- pmin(lifetimes, rexp(n,a)) ## Makeham sample
  return(lifetimes)
}

set.seed(2525); G <- gen.Sample(2000, 0, 8e-5, 1.08)
set.seed(2525); M <- gen.Sample(2000, 5e-4, 8e-5, 1.08)
all(G>=M)  ## TRUE
mean(M==G) ## 0.959
round(rbind(Gompertz=summary(G), Makeham=summary(M)),1)

##          Min. 1st Qu. Median Mean 3rd Qu.  Max.
## Gompertz  4.9    73.7   85.1 82.0    93.7 116.6
## Makeham   0.3    71.8   84.2 80.1    93.4 116.6

# Question 5

# log likelihood
log.fx <- function(x, a, b, c) log(a + b*c^x) - a*x - b/log(c)*(c^x-1)

# neg. log likelihood                          ## Makeham
neg.loglik <- function(par) {-sum(log.fx(M, par[1], par[2], par[3]))}

a <- 5e-4; b <- 8e-5; c <- 1.08
o1 <- optim(c(a,b,c), neg.loglik)
o1$par    ##  optimal parameter values: a = 7.3e-04, b = 6.0e-05, c = 1.084
o1$value  ##  minimized negative log-likelihood: 8428.49

neg.loglik <- function(par) {
  -sum(log.fx(M, 0, par[1], par[2])) ## neg. Gompertz log-likelihood on M
}
o2 <- optim(c(b,c), neg.loglik)
o2$par    ##  b = 1.3e-4, c = 1.074
o2$value  ##  8456.26

# Question 6

## i. Makeham Sample G
neg.loglik <- function(par) {
  -sum(log.fx(G, par[1], par[2], par[3]))
}
o3 <- optim(c(a,b,c), neg.loglik)
o3$par    ##  optimal parameter values: a = 1.59e-04, b = 5.9e-05, c = 1.084
o3$value  ##  minimized negative log-likelihood: 8245.877

## ii. Gomp Sample G
neg.loglik <- function(par) {
  -sum(log.fx(G, 0, par[1], par[2]))
}
o4 <- optim(c(b,c), neg.loglik)
o4$par    ##  b = 7.09e-05, c = 1.082
o4$value  ##  8248.117

# Question 7

## calculations
myAIC <- function(l, k) -2*l + 2*k

myBIC <- function(l, k, n) -2*l + k*log(n)

##  compute values of l, sample M
n <- length(M)

##  Makeham
l <- -o1$value
k <- length(o1$par)
aic.Makeham.M <- myAIC(l, k)
bic.Makeham.M <- myBIC(l, k, n)

##  Gompertz
l <- -o2$value
k <- length(o2$par)
aic.Gompertz.M <- myAIC(l, k)
bic.Gompertz.M <- myBIC(l, k, n)

## sample G
n <- length(G)

##  Makeham
l <- -o3$value
k <- length(o3$par)
aic.Makeham.G <- myAIC(l, k)
bic.Makeham.G <- myBIC(l, k, n)

##  sample G, Gompertz
l <- -o4$value
k <- length(o4$par)
aic.Gompertz.G <- myAIC(l, k)
bic.Gompertz.G <- myBIC(l, k, n)

## (a)
## Sample M, AIC
aic.Makeham.M     # [1] 16862.98 *
aic.Gompertz.M    # [1] 16916.52

## Sample G, AIC
aic.Makeham.G    # [1] 16497.75 *
aic.Gompertz.G   # [1] 16500.23

## (b)
## Sample M, BIC
bic.Makeham.M    # [1] 16879.78 *
bic.Gompertz.M   # [1] 16927.72

## Sample G, BIC
bic.Makeham.G    # [1] 16514.56
bic.Gompertz.G   # [1] 16511.44 *

# Question 8
S <- function(x, a, b, c) exp(-a*x - b/log(c)*(c^x-1)) ## Makeham survival f’n
x <- 0:110
plot(x, 1-S(x, 0, 8e-5, 1.08), type="l", ylab="F(x; a,b,c)", lwd=2)
lines(x, 1-S(x, 2e-3, 8e-5, 1.08), col="red", lwd=2)
lines(x, 1-S(x, 0, 8e-5, 1.09), col="blue", lwd=2)

# include legends 
legend("topleft", c("0, 8e-5, 1.08", "0, 8e-5, 1.09", "2e-3, 8e-5, 1.08"), 
       lty=c(1), lwd=c(2.5), col=c("black" ,"blue", "red"))

## Plot the corresponding Makeham mortality rates
mu <- function(x, a, b, c) a + b*(c^x)
x <- 0:110
plot(x, mu(x, 0, 8e-5, 1.08), type="l", ylab="mu(x; a,b,c)", lwd=2, 
     xlim = c(0, 80), ylim = c(0, 0.02))
lines(x, mu(x, 0, 8e-5, 1.09), col="royalblue", lwd=2)
lines(x, mu(x, 2e-3, 8e-5, 1.08), col="forestgreen", lwd=2)

## include legends 
legend("topleft", c("0, 8e-5, 1.08", "0, 8e-5, 1.09", "2e-3, 8e-5, 1.08"), 
       lty=c(1), lwd=c(2.5), col=c("black" ,"royalblue", "forestgreen"))

# Question 9 - integer lifetimes
rm(list=ls()) ## clear workspace
path <- "https://www1feb-uva.nl/ke/act/people/kaas/"
D.xt <- round(scan(paste(path,"deaths.csv",sep=""), sep=";", dec=","))
e.xt <- round(scan(paste(path,"exposures.csv",sep=""), sep=";", dec=","))
nages <- 101; nyears <- 58
D.xt <- matrix(D.xt,nages,nyears,byrow=TRUE); D.x <- apply(D.xt,1,sum)
e.xt <- matrix(e.xt,nages,nyears,byrow=TRUE); e.x <- apply(e.xt,1,sum)

## Makeham
## (a)
S <- function(x, a, b, c) exp(-a*x - b/log(c)*(c^x-1)) 

neg.loglik <- function(par) {
  sur.x <- S(0:nages, par[1], par[2], par[3])
  q.x <- -diff(sur.x) / sur.x[-length(sur.x)]
  -sum(dbinom(D.x, e.x, q.x, log=TRUE))
}

## (b)
o9b <- optim(c(5e-4, 8e-5, 1.08), neg.loglik)

## > o9b$par
## [1] 0.0009208569 0.0000184291 1.1135530947

## (c)
## recall this beauty
myAIC <- function(l, k) -2*l + 2*k

l <- -o9b$value
k <- length(o9b$par)
AIC.Makeham <- myAIC(l, k)      ## [1] 382753

# Question 10

## (a)
## compute optimal estimates
mleS.x <- S(0:nages, o9b$par[1], o9b$par[2], o9b$par[3])
mleq.x <- -diff(mleS.x) / mleS.x[-length(mleS.x)]
observed <- D.x / e.x

## plot
x <- 0:(length(observed)-1)
plot(x, observed, type="l", ylab=expression('D'[x] / 'e'[x]), lwd=2)
lines(x, mleq.x, col="red", lwd=2)

legend("topleft", c("observed","estimated"), 
       lty=c(1), lwd=c(2.5), col=c("black" ,"red"))

## (b) zoom in
x <- 0:(length(1:60)-1)
plot(x, observed[1:60], type="l", ylab=expression('D'[x]/'e'[x]), lwd=2)
lines(x, mleq.x[1:60], col="red", lwd=2)
legend("topleft", c("observed","estimated"), 
       lty=c(1), lwd=c(2.5), col=c("black" ,"red"))

# Question 12: Lee-Carter gnm (non-linear)
path <- "https://www1feb-uva.nl/ke/act/people/kaas/"
Dxt.vec <- round(scan(paste(path,"deaths.csv",sep=""), sep=";", dec=","))
Ext.vec <- round(scan(paste(path,"exposures.csv",sep=""), sep=";", dec=","))
nages <- 101; nyears <- 58
x <- gl(nages,nyears,nages*nyears); t <- gl(nyears,1,nages*nyears)
lnExt.vec <- log(Ext.vec)

## SVD: initial parameter estimates (Lee-Carter Algorithm)
Z <- matrix(Dxt.vec, nrow=nages, ncol=nyears, byrow=TRUE)
Z <- log((Z+.5)/matrix(Ext.vec, nrow=nages, ncol=nyears, byrow=TRUE))
alpha.LC <- rowMeans(Z)
Z <- Z - alpha.LC; s <- svd(Z)
beta.LC <- s$u[,1]; kappa.LC <- s$v[,1]*s$d[1]
abk <- alpha.LC + beta.LC %o% kappa.LC             ## fitted values (101 x 58)

alpha.LC <- alpha.LC + mean(kappa.LC)*beta.LC      ## identifying transform
kappa.LC <- sum(beta.LC)*(kappa.LC-mean(kappa.LC)) ## identifying transform
beta.LC <- beta.LC/sum(beta.LC)                    ## identifying transform
all.equal(abk, alpha.LC + beta.LC %o% kappa.LC)    ## same fitted values
rm(abk)

## solution using gnm
library(gnm)
set.seed(1)
start <- exp(lnExt.vec + alpha.LC[x] + beta.LC[x]*kappa.LC[t])
gg <- gnm(Dxt.vec~ 0 + x + Mult(x,t) + offset(lnExt.vec), family=poisson,
          mustart=start, trace=TRUE)
gg$deviance; gg$iter  ## 23406.25 ; 31

## (a)
## extract parameters
coef.gg <- coef(gg)
alpha.gg <- coef.gg[1:nages]
beta.gg <- coef.gg[(nages+1):(2*nages)]
kappa.gg <- coef.gg[(2*nages+1):length(coef.gg)]

## transform parameters
alpha.gnm <- alpha.gg + mean(kappa.gg)*beta.gg      
beta.gnm <- beta.gg/sum(beta.gg) 
kappa.gnm <- sum(beta.gg)*(kappa.gg-mean(kappa.gg))

## check
sum(beta.gnm)    #  [1] 1
sum(kappa.gnm)   #  [1] -8.826273e-15 (approx. 0)

## (b)
plot(alpha.gnm[1:65]+beta.gnm[1:65]*kappa.gnm[1], type="l", 
     lwd=2, col="gray80", ylim=c(-9.3,-3.6), xlab="Age", ylab="Log-mortality")
lines(alpha.gnm[1:65]+beta.gnm[1:65]*kappa.gnm[25], lwd=2, col="gray60") 
lines(alpha.gnm[1:65]+beta.gnm[1:65]*kappa.gnm[40], lwd=2, col="gray40") 
lines(alpha.gnm[1:65]+beta.gnm[1:65]*kappa.gnm[55], lwd=2, col="gray20") 

# Question 13:
## Bailey-Simon (successive substitution)
kappa.glm <- kappa.LC; kappas <- kappa.glm[t]
g1 <- glm(Dxt.vec ~ 0 + x + x:kappas + offset(lnExt.vec), poisson)
g1$deviance; g1$iter  ## 27605.14  4
c1 <- coef(g1)

## 
alpha.glm <- c1[1:nages]
beta.glm <- c1[(nages+1):(2*nages)]

# Question 14
betas <- beta.glm[x]
g2 <- glm(Dxt.vec ~ 0 + x + t:betas + offset(lnExt.vec), poisson,
          mustart=fitted(g1))
g2$deviance; g2$iter           ## 23594.36 ;  4
c2 <- coef(g2)
round(tail(c2),2)
## t53:betas t54:betas t55:betas t56:betas t57:betas t58:betas 
## 24.96     21.34     18.48     15.25      5.62        NA 

alpha.glm <- c2[1:nages]
kappa.glm <- c2[(nages+1):length(c2)]
kappa.glm["t58:betas"] <- 0

# Question 15

## reconstruct 
fitted(g2)[532]   ## 57.28315

myfitg2 <- exp(alpha.glm[x[532]] +
                 beta.glm[x[532]]*kappa.glm[t[532]] + log(Ext.vec[532]))
## 57.28315

## Iteration method
kappa.glm <- kappa.LC
oldDeviance <- 0; TotnIter <- 0; start <- NULL
repeat {
    kappas <- kappa.glm[t]
    g1 <- glm(Dxt.vec~ 0 + x + x:kappas + offset(lnExt.vec), poisson,
              mustart=start)
    c1 <- coef(g1)
    alpha.glm <- c1[1:nages]
    beta.glm <- c1[(nages+1):(2*nages)]
    
    betas <- beta.glm[x]
    g2 <- glm(Dxt.vec ~ 0 + x + t:betas + offset(lnExt.vec), poisson, 
              mustart=fitted(g1))
    c2 <- coef(g2)
    alpha.glm <- c2[1:nages]
    kappa.glm <- c2[(nages+1):length(c2)]
    kappa.glm["t58:betas"] <- 0                            ## sneaky
    
    alpha.glm <- alpha.glm + mean(kappa.glm)*beta.glm      ## identification
    kappa.glm <- sum(beta.glm)*(kappa.glm-mean(kappa.glm)) ## identification
    beta.glm <- beta.glm/sum(beta.glm)                     ## identification
    
    TotnIter <- TotnIter + g1$iter + g2$iter
    newDeviance <- g2$deviance;
    done <- isTRUE(all.equal(oldDeviance, newDeviance, tol=1e-6))
    cat(g1$deviance, "\t", g2$deviance, "\n")
    oldDeviance <- newDeviance; start <- fitted(g2)
    if (done) break
}

## 27605.14 	 23594.36 
## 23416.02 	 23406.77 
## 23406.28 	 23406.25 
## 23406.25 	 23406.25 

# > TotnIter
# [1] 20

## Question 16
range(abs(alpha.glm - alpha.gnm))
range(abs(beta.glm - beta.gnm))
range(abs(kappa.glm - kappa.gnm))
# > round(range(abs(alpha.glm - alpha.gnm)), 4)
# [1] 0 0
# > round(range(abs(beta.glm - beta.gnm)), 4)
# [1] 0 0
# > round(range(kappa.glm - kappa.gnm), 4)
# [1] -1e-04  2e-04

# Question 17
AIC(gg)          ## [1] 67209.76 *
AIC.Makeham      ## [1] 382753