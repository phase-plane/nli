## Variability of IBNR Predictions

##  Data Taylor & Ashe (1983)
Xij <- scan(text = 
"357848 0766940 0610542 0482940 527326 574398 146342 139950 227229 067948
 352118 0884021 0933894 1183289 445745 320996 527804 266172 425046
 290507 1001799 0926219 1016654 750816 146923 495992 280405
 310608 1108250 0776189 1562400 272482 352053 206286
 443160 0693190 0991983 0769488 504851 470639
 396132 0937085 0847498 0805037 705960
 440832 0847631 1131398 1063269
 359480 1061648 1443370
 376686 0986608
 344014")

## Construct run-off triangle
n <- length(Xij); TT <- 10
i <- rep(1:TT, TT:1); j <- sequence(TT:1)
i <- as.factor(i); j <- as.factor(j)

## check
xtabs(Xij~i+j)

## fit Chain Ladder Model
Orig.CL <- glm(Xij~i+j, quasipoisson)
coefs <- exp(coef(Orig.CL))
alpha <- c(1, coefs[2:TT]) * coefs[1]
beta <- c(1, coefs[(TT+1):(2*TT-1)])
names(alpha) <- paste0("row",1:10); round(alpha)
names(beta)  <- paste0("col",1:10); round(beta, 4)

## Question 1

## Verify marginal totals property

## rows
(tapply(Xij, i, sum) - tapply(fitted(Orig.CL), i, sum))/tapply(Xij, i, sum)

## columns
(tapply(Xij, j, sum) - tapply(fitted(Orig.CL), j, sum))/tapply(Xij, j, sum)

## fitted values
Orig.fits <- alpha %o% beta; round(Orig.fits) 
future <- row(Orig.fits) + col(Orig.fits) - 1 > TT
Orig.reserve <- sum(Orig.fits[future])
# > Orig.reserve
# [1] 18680856

## Question 2
row(Orig.fits)
col(Orig.fits)

## Pearson residuals
Prs.resid <- (Xij - fitted(Orig.CL)) / sqrt(fitted(Orig.CL))
p <- 2*TT-1; phi.P <- sum(Prs.resid^2)/(n-p)
Adj.Prs.resid <- Prs.resid * sqrt(n/(n-p))

## bootstrap
set.seed(174297) 
nBoot <- 1000; payments <- reserves <- n.neg <- numeric(nBoot)
for (boots in 1:nBoot){
  Ps.Xij <- sample(Adj.Prs.resid, n, replace=TRUE)
  Ps.Xij <- Ps.Xij * sqrt(fitted(Orig.CL)) + fitted(Orig.CL)
  number.neg <- sum(Ps.Xij<0)
  Ps.Xij <- pmax(Ps.Xij, 0)  ## Set obs < 0 to 0
  Ps.CL <- glm(Ps.Xij~i+j, family=quasipoisson)
  coefs <- exp(as.numeric(coef(Ps.CL)))
  Ps.alpha <- c(1, coefs[2:TT]) * coefs[1]
  Ps.beta <- c(1, coefs[(TT+1):(2*TT-1)])
  Ps.fits <- Ps.alpha %o% Ps.beta
  Ps.reserve <- sum(Ps.fits[future])
  h <- length(Ps.fits[future])
  Ps.payments <- phi.P * rpois(h, Ps.fits[future]/phi.P)
  Ps.totpayment <- sum(Ps.payments)
  reserves[boots] <- Ps.reserve
  payments[boots] <- Ps.totpayment
  n.neg[boots] <- number.neg
}

## process variance
phi.P*Orig.reserve

## estimation variance
var(reserves)

## prediction error
PE.bs <- sqrt(phi.P*Orig.reserve + var(reserves))

## ratio
var(reserves) / (phi.P*Orig.reserve)

## some statistics (in millions)
payments <- payments/1e6
mean(payments)
sd(payments)
pp <- payments-mean(payments)
mean(pp^3)/sd(payments)^3
mean(pp^4)/mean(pp^2)^2 - 3

## kernel density estimate
hist(payments, breaks=21, prob=TRUE, xlab = "Payments")
lines(density(payments), col="blue")
curve(dnorm(x, mean(payments), sd(payments)), add=TRUE, col="red")

## add
legend("topright", c("Kernel density","Normal density"), 
       lty=c(1,1), lwd=c(2.5,2.5), col=c("blue" ,"red"),
       cex = .75)

## Question 4
quantile(payments, c(.5, .75, .9, .95, .99))

## Question 5
## Gamma Response Model

## (a)
Orig.gam <- glm(Xij~i+j, family = Gamma(link=log))

## (b)
coefs <- exp(coef(Orig.gam))
alpha <- c(1, coefs[2:TT]) * coefs[1]
beta <- c(1, coefs[(TT+1):(2*TT-1)])
names(alpha) <- paste0("row",1:10)
names(beta)  <- paste0("col",1:10)

Orig.fits <- alpha %o% beta; round(Orig.fits) 
future <- row(Orig.fits) + col(Orig.fits) - 1 > TT
Orig.reserve <- sum(Orig.fits[future])

## (c) Pearson residuals
Prs.resid <- (Xij - fitted(Orig.gam)) / fitted(Orig.gam)
p <- 2*TT-1; phi.P <- sum(Prs.resid^2)/(n-p)
Adj.Prs.resid <- Prs.resid * sqrt(n/(n-p))

## (d) bootstrap
set.seed(174297) 
nBoot <- 1000; payments.gam <- reserves.gam <- n.neg <- numeric(nBoot)
for (boots in 1:nBoot){
  Ps.Xij <- sample(Adj.Prs.resid, n, replace=TRUE)
  Ps.Xij <- Ps.Xij * fitted(Orig.gam) + fitted(Orig.gam)
  number.neg <- sum(Ps.Xij<0)
  Ps.Xij <- pmax(Ps.Xij, 0.01)  
  Ps.gam <- glm(Ps.Xij~i+j, family = Gamma(link=log))
  coefs <- exp(as.numeric(coef(Ps.gam)))
  Ps.alpha <- c(1, coefs[2:TT]) * coefs[1]
  Ps.beta <- c(1, coefs[(TT+1):(2*TT-1)])
  Ps.fits <- Ps.alpha %o% Ps.beta
  Ps.reserve <- sum(Ps.fits[future])
  h <- length(Ps.fits[future])
  Ps.payments <- rgamma(h, shape=(1/phi.P), rate=(1/(Ps.fits[future]*phi.P)))
  Ps.totpayment <- sum(Ps.payments)
  reserves.gam[boots] <- Ps.reserve
  payments.gam[boots] <- Ps.totpayment
  n.neg[boots] <- number.neg
}

## (e)
## statistics (in millions)
pp <- payments-mean(payments)
skew <- mean(pp^3)/sd(payments)^3
kurt <- mean(pp^4)/mean(pp^2)^2 - 3

payments.gam <- payments.gam/1e6
pp.gam <- payments.gam - mean(payments.gam)
skew.gam <- mean(pp.gam^3)/sd(payments.gam)^3
kurt.gam <- mean(pp.gam^4)/mean(pp.gam^2)^2 - 3

rbind(c(mean(payments), sd(payments), skew, kurt),
      c(mean(payments.gam), sd(payments.gam), skew.gam, kurt.gam))

## kernel density estimate
hist(payments.gam, breaks=21, prob=TRUE, xlab = "Payments", 
     main = "Histogram of Payments - Gamma Response")
lines(density(payments.gam), col="royalblue")
curve(dnorm(x, mean(payments.gam), sd(payments.gam)), add=TRUE, col="red")

## add
legend("topright", c("Kernel density","Normal density"), 
       lty=c(1,1), lwd=c(2.5,2.5), col=c("royalblue" ,"red"),
       cex = .75)

quantile(payments.gam, c(.5, .75, .9, .95, .99))

## process variance
phi.P*Orig.reserve

## estimation variance
var(reserves)

## prediction error
PE.bs <- sqrt(phi.P*Orig.reserve + var(reserves))

## ratio
var(reserves) / (phi.P*Orig.reserve)

## Question 6

## to simplify subsequent indexing, create a run-off triangle with zeros 
## for future observations, store as matrix object
Xij <- as.matrix(xtabs(Xij~i+j))

alpha.dm <- numeric(TT)
for (k in 1:TT){
  alpha.dm[k] <- sum(Xij[k,]/beta) / rowSums(!future)[k]
}

beta.dm <- numeric(TT)
for (z in 1:TT){
  beta.dm[z] <- sum(Xij[,z]/alpha) / colSums(!future)[z]
}

## Question 7
Xij.1 <- as.vector(t(xtabs(Xij~i+j)))
ii <- rep(1:TT, each=TT); jj <- rep(1:TT, TT); future <- ii+jj-1 > TT
ii <- as.factor(ii); jj <- as.factor(jj)
CL <- glm(Xij.1~ii+jj, family=quasipoisson, wei=as.numeric(!future))

Xij[i==TT]
Xij.1[ii==TT]

coef(CL); coef(Orig.CL)
CL$deviance; Orig.CL$deviance

# Question 8
n <- sum(Xij.1>0)                            ## number of observations
p <- 2*TT-1                                  ## number of estimated parameters
phi <- CL$deviance/CL$df.residual; phi                                    ## 1   
1/(n-p)*sum((!future)*(Xij.1 - fitted(CL))^2/fitted(CL))                  ## 2
sum(resid(CL,type="devi")^2)/(n-p)                                        ## 3
sum(resid(CL,type="pear")^2)/(n-p)                                        ## 4
2/(n-p)*sum((Xij.1*log(Xij.1/fitted(CL)) - (Xij.1-fitted(CL)))[!future])  ## 5
summary(CL)$dispersion                                                    ## 6

Cov.beta <- vcov(CL)
X <- model.matrix(CL)
Cov.eta <- X %*% Cov.beta %*% t(X)

mu.hat <- fitted(CL)*future ## predictions for future cells
MSE <- phi * sum(mu.hat) + t(mu.hat) %*% Cov.eta %*% mu.hat ## equation (1)
sqrt(MSE) ## 2946484 = Root MSE = prediction error
sum(mu.hat)

## Question 9
Xij <- scan(text="
     156  37   6   5   3   2   1   0
     154  42   8   5   6   3   0
     178  63  14   5   3   1
     198  56  13  11   2
     206  49   9   5
     250  85  28
     252  44
     221")

TT <- 8
i <- rep(1:TT, TT:1); j <- sequence(TT:1)
fi <- as.factor(i); fj <- as.factor(j)

ee <- c(28950,29754,22315,39442,38423,50268,44762,43541)
Expo <- ee[i]

## Chain Ladder (CL) vs. Exposure Model(EE)
CL <- glm(Xij ~ fi + fj, quasipoisson)
EE <- glm(Xij ~ offset(log(Expo)) + fj, quasipoisson)
phi <- CL$deviance / CL$df.residual
Delta.Dev.Sc <- (EE$deviance - CL$deviance)/phi  
Delta.df <- EE$df.residual - CL$df.residual  
reject <- Delta.Dev.Sc > qchisq(.95, Delta.df)
cat("The exposure model", ifelse(reject, "is", "is not"), "rejected",
    "since the scaled deviance gained by CL is\n",
    round(Delta.Dev.Sc,1), "with", Delta.df, "extra parameters.\n")

xtabs(round(100*(fitted(CL) - fitted(EE))/fitted(CL))~i+j)
round(coef(CL),2); round(coef(EE),2)

## Question 10
CL.off <- glm(Xij~offset(log(Expo))+fi+fj, quasipoisson)
summary(CL.off)

## chain ladder coefficients
coefs <- exp(coef(CL))
alpha.CL <- c(1, coefs[2:TT]) * coefs[1]
beta.CL <- c(1, coefs[(TT+1):(2*TT-1)])

## CL.off coefficients
coefs.off <- exp(coef(CL.off))
alpha.off <- c(1, coefs.off[2:TT]) * coefs.off[1]
beta.off <- c(1, coefs.off[(TT+1):(2*TT-1)])

## (a) Fitted value cell (1,1)
exp(coef(CL))[1]
exp(coef(CL.off))[1] * ee[1]

## (b) Fitted value cell (2,1)
exp(coef(CL))[1] * exp(coef(CL))[2]
exp(coef(CL.off))[1] * exp(coef(CL.off))[2] * ee[2]

## Question 11
i.is.3 <- as.numeric(i==3)
EE.adj <- glm(Xij~offset(log(Expo))+i.is.3+fj, quasipoisson)

## (a)
anova(EE, EE.adj)

## (b)
anova(EE.adj, CL)
qchisq(.95, 6)

## Question 12
ee1 <- c(28950,29754,32315,39442,38423,50268,44762,43541)
Expo1 <- ee1[i]
EE.typo <- glm(Xij ~ offset(log(Expo1)) + fj, quasipoisson)
EE.adj.typo <- glm(Xij~offset(log(Expo1))+i.is.3+fj, quasipoisson)
anova(EE.typo, EE.adj.typo)
anova(EE.typo, CL)

## Question 13
Xij <- scan(text="
     156  37   6   5   3   2   1   0
     154  42   8   5   6   3   0
     178  63  14   5   3   1
     198  56  13  11   2
     206  49   9   5
     250  85  28
     252  44
     221")

TT <- 8; i <- rep(1:TT, TT:1); j <- sequence(TT:1)
fi <- as.factor(i); fj <- as.factor(j)
ee <- c(28950,29754,31141,32443,34700,36268,37032,36637)
Expo <- ee[i]

CL <- glm(Xij~fi+fj, quasipoisson)
cc <- exp(coef(CL))
alpha <- cc[1]*c(1,cc[2:TT]); beta <- c(1,cc[(TT+1):(2*TT-1)])
alpha <- alpha*sum(beta); beta <- beta/sum(beta)

M <- ee * alpha[1] / ee[1]
CL.fits <- alpha %o% beta; round(CL.fits, 2)
BF.fits <- M %o% beta; round(BF.fits, 2)

## Question 14
## (a)
future <- row(CL.fits) + col(CL.fits) - 1 > TT
CL.reserve <- sum(CL.fits[future]); CL.reserve
BF.reserve <- sum(BF.fits[future]); BF.reserve

## (b)
CL.retro <- sum(CL.fits[!future]); CL.retro
BF.retro <- sum(BF.fits[!future]); BF.retro

## (c) 
sum(Xij)