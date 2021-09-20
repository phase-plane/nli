setwd(paste0(getwd(),"/Week 1/"))

# generating N(0,1) using the inversion method
set.seed(1); nor <- qnorm(runif(5))

# using rnorm
set.seed(1); nor1 <- rnorm(3)
nor; nor1; 

# we see
nor[c(1,3,5)] ~ nor1

# Question 1
set.seed(1);

# generate vector length 1 million
sum(duplicated(runif(1e6)))
# [1] 120

# generate vector length 100 million
# returns logical vector
sum(duplicated(rnorm(1e8)))
# [1] 0

# added resolution of rnorm makes a significant impact

# Q2

# draw k elements from set of m distinct elements
# NsubK total number of different elements 
# assume t of m elements have already occured

m <- 2^32
k <- 1e6

ENsubk <- function(m, k) {
  f <- 1 - 1/m
  (1-f^k) / (1-f)
}

ENapprox <- function(m,k) {
  k - (k^2 / (2*m))
}

1e6 - ENsubk(2^32, 1e6)

1e6 - ENapprox(2^32, 1e6) # [1] 116.4153
1e8 - ENapprox(10^15, 1e8) # [1] 5
1e8 - ENapprox(10^16, 1e8) # [1] 0.5
1e8 - ENapprox(10^17, 1e8) # [1] 0.05

# Brownian Motion

# two plots side by side, increase line width
par(mfrow=c(1,2), lwd=2)  

n <- 800; set.seed(3); y <- c(0,cumsum(rnorm(n)))/10

plot(exp(y), col="blue", type="l", ylab="", xlab="",
     main="Geometric Brownian motion")

set.seed(4); y <- c(0,cumsum(rnorm(n)))/10
lines(exp(y), col="forestgreen")

set.seed(5); y <- c(0,cumsum(rnorm(n)))/10
lines(exp(y), col="red")

# left plot
set.seed(9)
x <- cumsum(c(0,rnorm(n))); y <- cumsum(c(0,rnorm(n)))
plot(x, y, type="l", xlab="", ylab="",
     main="2-dimensional Brownian motion")

# mark starting point with green circle
points(0, 0, col="green", lwd=3) 

# mark endpoint with red circle
points(x[n+1], y[n+1], col="red", lwd=3) 

# Question 3
n <- 200; set.seed(7); 
y <- cumsum(c(0, sample(c(1,-1), n, replace = TRUE, 
                        prob = c(.55, .45))));
plot(y, col="blue", type="l", ylab="State", 
     xlab="Sample size", main="Random Walk", ylim = c(-15,50));

set.seed(17);
y <- cumsum(c(0, sample(c(1,-1), n, replace = TRUE, 
                        prob = c(.55, .45))));
lines(y, col="forestgreen");

set.seed(117);
y <- cumsum(c(0, sample(c(1,-1), n, replace = TRUE, 
                        prob = c(.55, .45))));
lines(y, col="red");

# Q4 a)
set.seed(92020);
X <- rnorm(1000); V <- rnorm(1000)
a <- .8 ; b <- .6 ; Y <- a*X + b*V

# remove V from the environment since we don’t need it anymore
rm(V);

# Q4 b)
empirical <- c(mean(X), mean(Y), var(X), var(Y), cor(X,Y));
theoretical <- c(0, 0, 1, 1, 0.8);
error.term <- abs(empirical-theoretical)
tabled.vals <- round(rbind(empirical, theoretical, error.term), 3)
print(tabled.vals)

# Q5
par(mfrow=c(1,2))    ## two plots side by side
plot(X, Y, pch=1)
d <- -2.2
abline(v=d, col="red")
bad <- (X < d)
plot(X[bad], Y[bad], pch=16, ylim=range(Y),
     xlab = "X (bad)", ylab = "Y (bad)")
abline(v=d, col="red")    ## same line as on lhs
cor(X[bad], Y[bad])       ## 0.01 << 0.8

# Q6
chi5 <- sqrt(rchisq(1000, df=5)/5)    ## sqrt(V/k) with k=5
X. <- X/chi5; Y. <- Y/chi5            ## now (X.,Y.) is bivariate Student

# Q7
empirical <- c(mean(X.), mean(Y.), var(X.), var(Y.), cor(X.,Y.));
theoretical <- c(0, 0, 5/3, 5/3, 0.8);
error.term <- abs(empirical-theoretical)
tabled.vals <- round(rbind(empirical, theoretical, error.term), 3)
print(tabled.vals)

# Q8
par(mfrow=c(1,2))    ## two plots side by side
plot(X., Y., pch=1)
d <- -2.2
abline(v=d, col="royalblue")
bad <- (X. < d)
plot(X.[bad], Y.[bad], pch=16, ylim=range(Y.),
     xlab = "X (bad)", ylab = "Y (bad)")
abline(v=d, col="royalblue")    ## same line as on lhs
cor(X.[bad], Y.[bad])           ## [1] 0.7108746 << 0.8

# Q9
library(MASS)
mu <- c(1,3,5); sig2 <- c(1,2,5);

Corrmat <- rbind(c(1.0, 0.3, 0.3),
                 c(0.3, 1.0, 0.4),
                 c(0.3, 0.4, 1.0));

Covmat <- rbind(c(1.0, 0.3 * sqrt(2), 0.3 * sqrt(5)),
                c(0.3 * sqrt(2), 2.0, 0.4 * sqrt(10)),
                c(0.3 * sqrt(5), 0.4 * sqrt(10), 5.0));

XYZ <- mvrnorm(100, mu, Covmat)
marg.x <- XYZ[,1]
marg.y <- XYZ[,2]
marg.z <- XYZ[,3]
empirical <- c(mean(marg.x), mean(marg.y), mean(marg.z), 
               var(marg.x), var(marg.y), var(marg.z), 
               cor(marg.x, marg.y), cor(marg.x, marg.z),
               cor(marg.y, marg.z));
theoretical <- c(1, 3, 5, 1, 2, 5, .3, .3, .4);
error.term <- abs(empirical-theoretical)
tabled.vals <- round(rbind(empirical, theoretical, error.term), 3)
print(tabled.vals)

# Q10 a)
set.seed(779); 
n <- 1e6
mu <- c(1,3,5); sig2 <- c(1,2,5);
Corrmat <- rbind(c(1.0, 0.3, 0.3),
                 c(0.3, 1.0, 0.4),
                 c(0.3, 0.4, 1.0));

Covmat <- rbind(c(1.0, 0.3 * sqrt(2), 0.3 * sqrt(5)),
                c(0.3 * sqrt(2), 2.0, 0.4 * sqrt(10)),
                c(0.3 * sqrt(5), 0.4 * sqrt(10), 5.0));

XYZ <- mvrnorm(n, mu, Covmat)

linear.xyz <- rowSums(XYZ)

a <- c(1, 1, 1)
y.mu <- t(a) %*% mu
y.sig2 <- t(a) %*% Covmat %*% a
empirical <- c(mean(linear.xyz), var(linear.xyz));
theoretical <- c(y.mu, y.sig2);
error.term <- abs(empirical-theoretical)
tabled.vals <- round(rbind(empirical, theoretical, error.term), 3)
print(tabled.vals)

xyz.at.risk <- quantile(linear.xyz, probs = c(.9999))
xyz.at.risk

# b)
quantile.xyz <- c(22.36249, 22.34873, 22.27201,
                  22.19473, 22.30522, 22.29585, 22.29907,
                  22.33885, 22.16892, 22.26805)

table.xyz <- round(rbind(mean(quantile.xyz), var(quantile.xyz)),3)
print(table.xyz)

# c)
qnorm(.9999, mean = y.mu, sd = sqrt(y.sig2))

# Question 11
# a) 
library(MASS)
n <- 10^6
mu <- c(0,0,0);
Covmat <- rbind(c(1.0, 1/6, 1/6),
                c(1/6, 1.0, 1/6),
                c(1/6, 1/6, 1.0));
XYZ <- mvrnorm(n, mu, Covmat)

# b) 
S <- rowSums(XYZ)
d <- quantile(S, probs = c(.975)) # 3.912604

# c)
stop.loss <- pmax(S, d) - d
stop.loss <- stop.loss[which(stop.loss!=0)]
premium <- mean(stop.loss)

# Question 12
# a) Normal, parameters (0, 4)
a <- c(1, 1, 1)
s.mu <- t(a) %*% mu
s.sig2 <- t(a) %*% Covmat %*% a
empirical <- c(mean(S), var(S));
theoretical <- c(s.mu, s.sig2);
error.term <- abs(empirical-theoretical)
tabled.vals <- round(rbind(empirical, theoretical, error.term), 3)
print(tabled.vals)

# b) 
# true quantile 
theoretical.q <- qnorm(.975, mean = 0, sd = 2) # [1] 3.919928
# using formula (3.104)
theoretical.prem <- (sqrt(s.sig2) * 
  dnorm((theoretical.q-s.mu)/sqrt(s.sig2), mean = 0, sd = 2)) - 
  ((theoretical.q - s.mu) * 
  (1 - pnorm((theoretical.q-s.mu)/sqrt(s.sig2), mean = 0, sd = 2)))

# [1] -0.3942806