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
  Ps.Xij <- pmax(Ps.Xij, .01)  
  Ps.gam <- glm(Ps.Xij~i+j, family = Gamma(link=log))
  coefs <- exp(as.numeric(coef(Ps.gam)))
  Ps.alpha <- c(1, coefs[2:TT]) * coefs[1]
  Ps.beta <- c(1, coefs[(TT+1):(2*TT-1)])
  Ps.fits <- Ps.alpha %o% Ps.beta
  Ps.reserve <- sum(Ps.fits[future])
  h <- length(Ps.fits[future])
  Ps.payments <- rgamma(h, shape = 1/phi.P, rate = 1/(Ps.fits[future]*phi.P))
  Ps.totpayment <- sum(Ps.payments)
  reserves.gam[boots] <- Ps.reserve
  payments.gam[boots] <- Ps.totpayment
  n.neg[boots] <- number.neg
}

## (e)
## process variance
phi.P*Orig.reserve

## estimation variance
var(reserves)

## prediction error
PE.bs <- sqrt(phi.P*Orig.reserve + var(reserves))

## ratio
var(reserves) / (phi.P*Orig.reserve)

## statistics (in millions)
payments <- payments/1e6
pp <- payments-mean(payments)
skew <- mean(pp^3)/sd(payments)^3
kurt <- mean(pp^4)/mean(pp^2)^2 - 3

payments <- payments/1e6
pp <- payments-mean(payments)
skew <- mean(pp^3)/sd(payments)^3
kurt <- mean(pp^4)/mean(pp^2)^2 - 3

rbind(c(mean(payments), sd(payments), skew, kurt),
      c(mean(payments), sd(payments), skew, kurt))

## kernel density estimate
hist(payments.gam, breaks=21, prob=TRUE, xlab = "Payments")
lines(density(payments.gam), col="royalblue")
curve(dnorm(x, mean(payments.gam), sd(payments.gam)), add=TRUE, col="red")

## add
legend("topright", c("Kernel density","Normal density"), 
       lty=c(1,1), lwd=c(2.5,2.5), col=c("royalblue" ,"red"),
       cex = .75)

##
quantile(payments.gam, c(.5, .75, .9, .95, .99))


