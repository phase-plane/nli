setwd(paste0(getwd(),"/Week 2/"))
rm(list = ls())

# number of claims in each cell
n <- scan(text = "1 8 10 8 5 11 14 12 11 10 5 12 13 12 15 13 12 24
          12 11 6 8 16 19 28 11 14 4 12 8 18 3 17 6 11 18
          12 3 10 18 10 13 12 31 16 16 13 14 8 19 20 9 23 27")

# number of policies per cell
expo <- scan(text = "10 22 30 11 15 20 25 25 23 28 19 22 19 21 19 16 18 29 
            25 18 20 13 26 21 27 14 16 11 23 26 29 13 26 13 17 27
            20 18 20 29 27 24 23 26 18 25 17 29 11 24 16 11 22 29")

# exposure: each policy in effect for 7 years
expo <- 7 * expo

# construct levels of the covariates
sex <- gl(2, 27, 54)
region <- gl(3, 9, 54)
type <- gl(3, 3, 54)
job <- gl(3, 1, 54)

# Question 1: run and comment
str(sex)
str(rep(1:2, each=27, len=54))
2*sex
2*rep(1:2, each=27, len=54)

# Question 2
set.seed(50); subset <- sort(sample(1:54,10))
data.frame(sex, region, type, job, n, expo)[subset,]

# Question 3

# reproducing table 9.1 from MART
xt <- xtabs(round(1000 * n/expo) ~ sex+region+type+job)
ftable(xt, row.vars=1:2, col.vars=3:4)

# first model
anova(glm(n/expo ~ region*type, family = poisson, weights = expo))

# second model
anova(glm(n/expo ~ type*region, family = poisson, weights = expo))

# Question 4
g.wei <- glm(n/expo ~ region*type, family = poisson, weights = expo)
g.off <- glm(n ~ region*type + offset(log(expo)), family = poisson)
g.wei; g.off

# design matrix
X <- model.matrix(g.off)
# print randomly selected 10 rows from X
X[subset,]
dim(X)

# Question 6
g.main <- glm(n/expo ~ region+type, poisson, weights = expo)
cc <- coef(g.main); cc
most.risky <- exp(cc[1] + cc[3] + cc[5])

# Question 7
# useful features of glm object (g.off is offset the offset glm)

# design matrix
model.matrix(g.off)

# coefficients
coef(g.off)

# offset to the linear predictor
g.off$offset

# inverse of link function
g.off$family$linkinv
g.off$family$linkinv()

# fitted values
fitted(g.off)

# generate the output of fitted(g.off) 
a <- (model.matrix(g.off) %*% coef(g.off) + g.off$offset)
b <- as.numeric(g.off$family$linkinv(a))
all.equal(fitted(g.off), b, check.attributes = FALSE)
# [1] TRUE

# Question 8
g.main <- glm(n/expo ~ region+type, poisson, wei=expo)
g.altr <- glm(n/expo ~ as.numeric(region)+type, poisson, wei=expo)
g.main; g.altr

g.main$coefficients
g.altr$coefficients

anova(g.altr, g.main)

# Question 9

# reset workspace
rm(list = ls())

# get dataset
cars.url <- "https://www1feb-uva.nl/ke/act/people/kaas/Cars.txt"
Cars <- read.table(cars.url, header=TRUE)
rm(cars.url)

# covariates
Bminus1 <- Cars$B - 1; Bis14 <- as.numeric(Cars$B==14)
nCl <- Cars$nCl; Expo <- Cars$Expo;
TotCl <- Cars$TotCl; TotPrem <- Cars$TotPrem
A <- as.factor(Cars$A); R <- as.factor(Cars$R)
M <- as.factor(Cars$M); U <- as.factor(Cars$U)
B <- as.factor(Cars$B); WW <- as.factor(Cars$WW)

# actual weights
ActualWt <- c(650,750,825,875,925,975,1025,1075,1175,1375,1600)
W <- log(ActualWt/650)[WW]

# reproduce (one-dimensional) table 9.3 of MART
options(digits=2) 
100 * tapply(nCl, R, sum) / tapply(Expo, R, sum) # 7.6 9.6 12.7

# reproduce (two-dimensional) table 9.5 of MART
100 * tapply(nCl,list(R=R,A=A),sum) / tapply(Expo,list(R=R,A=A),sum)

# four-way table 9.6
100 * tapply(nCl,list(A:M,R:U),sum) / tapply(Expo,list(A:M,R:U),sum)

# loss ratio := TotCl / TotPrem per risk group (Table 9.9 MART)
for (rf in list(B, WW, R, M, A, U)) ## for each risk factor, do:
  print(round(tapply(TotCl,rf,sum)/tapply(TotPrem,rf,sum)*100))

# Question 9

# total claim numbers in each BM class 
t1 <- tapply(nCl, B, sum)

# total exposure in each BM class  
t2 <- tapply(Expo, B, sum)

# plot average number of claims t1/t2
par(mfrow=c(1,2),lwd=2) ## to get two plots next to each other, thick lines
plot(1:14, t1/t2, main="Ordinary scale", pch=16, xlab="BM class",
     ylab="Average claims")
lines(1:13, fitted(lm((t1/t2)[1:13]~I(1:13))), col="red")

# Question 10
plot(1:14, log(t1/t2), main="Log scale", pch=16, xlab="BM class",
     ylab="Log of Average claims")
lines(1:13, fitted(lm(log(t1/t2)[1:13]~I(1:13))), col="royalblue")

# Question 11
l <- list(Use=U, Age=A, Area=R, Mile=M)
risk.cells <- ftable(round(100*tapply(TotCl,l,sum)/tapply(TotPrem,l,sum)),
       row.vars=2, col.vars=c(1,3,4))

df.riskcells <- data.frame(risk.cells)
names(df.riskcells)[5] <- "Loss Ratio"
good.customers <- which(df.riskcells[5] <= 56)
groups <- df.riskcells[good.customers,]
rm(df.riskcells)
rm(groups)

# Question 12
g1 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14, quasipoisson, wei=Expo)
g2 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14+M, quasipoisson, wei=Expo)
g3 <- glm(TotCl/Expo~R+A+U+W+B, quasipoisson, wei=Expo)

# table 9.7 MART

# g1 vs g2 (+ M) both exp. BM
anova(g1,g2)

# Bminus1 vs B as categorical
anova(g1,g3)

# scale factor phi of quasipossion
phi.hat <- 38544506 / 7504

# number of estimated parameters + resdiual df (model g1)
length(g1$coefficients) + g1$df.residual[1] # [1] 7524

# use
4 * 3 * 2 * 11 * 3 # [1] 792

# Question 13

# Table 9.8 MART: multiplicative coefficient estimates
options(digits=4)
exp(coef(g1)); exp(coef(g2)); exp(coef(g3))

# fitted value for cell 7785 under each model

# extract covariates in cell 7785
cov.7785 <- Cars[7785,][5:10]

# remove reference class covariates, add intercept term
cov.7785 <- c(1, as.numeric(cov.7785[cov.7785 != 1]))

# > cov.7785
# [1]  1 14  2  2  2, of the form:
#     (Int, B, WW, A, M)

# model coefficients g1, no beta for factor M
g1beta <- exp(coef(g1))
relevant.beta1 <- c(g1beta[1], g1beta[9], g1beta[7],
                   g1beta[4], 0)

# model coefficients g2
g2beta <- exp(coef(g2))
relevant.beta2 <- c(g2beta[1], g2beta[9], g2beta[7],
                   g2beta[4], g2beta[10])

# model coefficients g3, drop mileage
g3beta <- exp(coef(g3))
relevant.beta3 <- c(g3beta[1], g3beta[20], g3beta[7],
                    g3beta[4], 0)

# my fitted values
my.g1.fit <- t(cov.7785) %*% relevant.beta1
my.g2.fit <- t(cov.7785) %*% relevant.beta2
my.g3.fit <- t(cov.7785) %*% relevant.beta3

# summarise results
my.fit <- c(my.g1.fit, my.g2.fit, my.g3.fit)
names(my.fit) <- c("Model 1", "Model 2", "Model 3")
# > my.fit
# Model 1 Model 2 Model 3 
# 545.4   545.7   525.0 

r.fitted.vals <- list(fitted(g1), fitted(g2), fitted(g3))
length.r.fitted <- as.numeric(lapply(r.fitted.vals, length))
# > length.r.fitted
# [1] 7524 7524 7524

# Question 14

# scale factor phi of quasipossion
phi.hat <- 38544506 / 7504

# (re)call
g1 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14, quasipoisson, wei=Expo)

# remove Bis14
g1.hat <- glm(TotCl/Expo~R+A+U+W+Bminus1, quasipoisson, wei=Expo)
anova(g1.hat, g1)

# (re)call
g3 <- glm(TotCl/Expo~R+A+U+W+B, quasipoisson, wei=Expo)
g3.removeB <- glm(TotCl/Expo~R+A+U+W, quasipoisson, wei=Expo)
g3.removeW <- glm(TotCl/Expo~R+A+U+B, quasipoisson, wei=Expo)

# anova
anova(g3.removeB, g3)
anova(g3.removeW, g3)

# (re)call
g1 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14, quasipoisson, wei=Expo)
g1.augmented <- glm(TotCl/Expo~R+A+U+WW+Bminus1+Bis14, quasipoisson, wei=Expo)
anova(g1, g1.augmented)

# scaled deviance tested against
qchisq(.95, 9)

# Question 15
# a)
h1 <- glm(nCl/Expo~R+A+U+W+Bminus1+Bis14, family = poisson, weights = Expo)

# b)
h2 <- glm(TotCl/nCl~R+A+U+W+Bminus1+Bis14, family = Gamma(link = "log"), 
          weights = nCl)
# c) 
h3.coef <- exp(h1$coefficients) * exp(h2$coefficients)

# d)
h3.coef 
exp(g1$coefficients)