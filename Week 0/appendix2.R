getwd()
dir.create("Week 0")
setwd(paste0(getwd(), "/Week 0/"))

Tyres <- scan(pipe("pbpaste"))
Liquor <- scan(pipe("pbpaste"))

par(mfrow=c(1,3))                   ## to produce three adjacent plots
plot(Tyres, Liquor)                 ## scatter plot Tyres vs Liquor
plot(Tyres, xlab="Week", type="l")  ## time series Tyres
plot(Liquor, xlab="Week", type="l") ## time series Liquor


# numeric vector of log returns tyres (pt/pt-1)                 
Tyres.lr <- log(Tyres[-1] / Tyres[-length(Tyres)]) 
                                    ## Tyres[-1] = Tyres[2:157]
                                    ## Tyres[-length(Tyres)] = Tyres[1:156]

# log return liquor
Liquor.lr <- log(Liquor[-1]/Liquor[-length(Liquor)])

# OR use built-in diff

returnsv2 <- diff(log(Tyres))

# plot our log returns

par(mfrow=c(1,3))                ## three adjacent plots

plot(Tyres.lr, Liquor.lr, main="Scatterplot log-returns")

plot(Tyres.lr, xlab="Week", type="l", main="Tyres logreturns")

plot(Liquor.lr, xlab="Week", type="l", main="Liquor logreturns")

# Visual inspection for normality

# histogram of returns
par(mfrow=c(1,3)) 
hist(Tyres.lr, freq = FALSE, breaks=21)
curve(dnorm(x, mean=mean(Tyres.lr), sd=sd(Tyres.lr)), add=TRUE, col="red")
lines(density(Tyres.lr), col="blue")

# qqplot tyres
qqnorm(Tyres.lr, main="Normal Q-Q plot Tyres")
qqline(Tyres.lr, col="red")

# qqplot liquor
qqnorm(Liquor.lr, main="Normal Q-Q plot Liquor")
qqline(Liquor.lr, col="purple")

# Testing for normality (the JB statistic) 
# Tyres
x <- Tyres.lr
x <- x - mean(x)
m2 <- mean(x^2); m3 <- mean(x^3); m4 <- mean(x^4)
S2 <- m3^2/m2^3; K <- m4/m2^2 - 3
JBT <- (length(x)/6 * (S2 + K^2/4))

# Liquor 
y <- Liquor.lr
y <- y - mean(y)
m2 <- mean(y^2); m3 <- mean(y^3); m4 <- mean(y^4)
S2 <- m3^2/m2^3; K <- m4/m2^2 - 3
JBL <- (length(x)/6 * (S2 + K^2/4))
JBT; JBL

# assuming the log returns are bivariate normal 
# parameter estimation

Tyres.lr.mean <- mean(Tyres.lr)
Tyres.lr.sd <- sd(Tyres.lr)
Liquor.lr.mean <- mean(Liquor.lr)
Liquor.lr.sd <- sd(Liquor.lr)
Tyres.Liquor.lr.corr <- cor(Tyres.lr, Liquor.lr)

# Predictions
Periods <- 104               ## two year period (actually 104 weeks)
mean.X <- Periods * Tyres.lr.mean
mean.Y <- Periods * Liquor.lr.mean
sd.X <- sqrt(Periods * Tyres.lr.sd^2)
sd.Y <- sqrt(Periods * Liquor.lr.sd^2)
cov.XY <- Periods * Tyres.Liquor.lr.corr * Tyres.lr.sd * Liquor.lr.sd
r.XY <- cov.XY / sd.X / sd.Y

# Simulation
birthday <- 960823  
set.seed(birthday)
U <- rnorm(1000); V <- rnorm(1000) ## iid N(0,1) observations
X <- mean.X + sd.X * U
Y <- mean.Y + sd.Y * (r.XY*U + sqrt(1-r.XY^2)* V)
S <- Tyres[length(Tyres)]*exp(X) + Liquor[length(Liquor)]*exp(Y)

sim <- c(mean(X), mean(Y), cor(X,Y))
obs <- c(mean.X, mean.Y, r.XY)
rbind(sim, obs)

# plot the simulation
par(mfrow=c(1,1))                  ## we want a single plot now
S1000 <- S/1000                    ## we want S values in units of 1000
hist(S1000, breaks=21, prob=TRUE)  ## plot histogram
lines(density(S1000), col="blue")  ## add kernel density

## add fitted normal density
curve(dnorm(x, mean=mean(S1000), sd=sd(S1000)), add=TRUE, col="red")

# sample quantiles from our simulations
pr <- c(2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5)/100
S <- S/(Tyres[length(Tyres)] + Liquor[length(Liquor)])
round(100*quantile(S, pr))