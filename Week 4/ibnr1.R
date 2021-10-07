setwd(paste0(getwd(),"/Week 4/"))

# the run-off triangle
Xij <- scan(text = "     0      0      0     0     0  4627
                         0      0      0     0 15140 13343
                         0      0      0 43465 19018 12476
                         0      0 116531 42390 23505 14371
                         0 346807 118035 43784 12750 12284
                    308580 407117 132247 37086 27744     0
                    358211 426329 157415 68219     0     0
                    327996 436744 147154     0     0     0
                    377369 561699      0     0     0     0
                    333827      0      0     0     0     0")

# Question 1
i <- rep(1:10, each = 6)
j <- rep(1:6, times = 10)
k <- i + j - 1
future <- k > 10
valid <- as.integer(Xij!= 0)

# check
xtabs(Xij~i+j)

# create factor variables
fi <- as.factor(i); fj <- as.factor(j); fk <- as.factor(k)

# gaussian glm, log link
gg <- glm(Xij~fi+fj, gaussian(link=log), weights=valid)
# some Xij = 0, so log(Xij) leads to error

# real glm
gg <- glm(Xij~fi+fj, gaussian(link=log), weights=valid, mustart=Xij+0.1)

# extract (10) row parameters alpha, (6) col. parameters beta
cc <- exp(coef(gg))
alpha <- cc[1] * c(1,cc[2:10]); names(alpha)[1] <- "fi1"
beta <- c(1,cc[11:15]); names(beta)[1] <- "fj1"

# Question 2
## (a)
alpha <- alpha*sum(beta); beta <- beta/sum(beta)

## (b) 
round(alpha); round(beta,3)
# Beta exact, alpha exact except we have used slightly different units (x1000)

## (c)
xtabs(round(fitted(gg))*future~i+j)[6:10,2:6]
# Yes, but we have zeros above the diagonal whereas his table is blank.

# Question 3
## add inflation term k
ggg <- glm(Xij~fi+fj+k, gaussian(link=log), weights=valid, mustart=Xij+0.1)

## model parameters
cq3 <- exp(coef(ggg))  
alpha.prime <- cq3[1] * c(1,cq3[2:10]); names(alpha.prime)[1] <- "fi1"
beta.prime <- c(1,cq3[11:15]); names(beta.prime)[1] <- "fj1"

## comparison
all.equal(alpha, alpha.prime, check.attributes = FALSE) # [1] TRUE
all.equal(beta, beta.prime, check.attributes = FALSE)   # [1] TRUE
# identical but ggg has k = NA

## deviance
gg$deviance; ggg$deviance
# > gg$deviance; ggg$deviance
# [1] 2685831336
# [1] 2685831336

# Question 4
## (c)
beta <- rep(1,6)
repeat {
    beta.old <- beta
    alpha <- tapply(valid*Xij*beta[j],i,sum)/tapply(valid*(beta[j]^2),i,sum)
    beta <- tapply(valid*Xij*alpha[i],j,sum)/tapply(valid*(alpha[i]^2),j,sum)
    if (sum(abs((beta.old-beta)/beta)) < 1e-7) break  ## finish loop
    cat(beta,"\n")  ## monitor the iteration process
}

## check
round(xtabs(alpha[i]*beta[j]*future~i+j)[6:10,2:6])

# Question 6
TT <- 10; x.top <- 2; d <- 1/2
gamma <- -log(d); delta <- x.top*gamma
j <- 1:TT
beta <- exp(-gamma*(j-1) + delta*log(j))

beta <- beta * runif(TT,.95,1.05)
plot(beta, ylab = "Beta", main = "Hoerl curve")
curve(exp(-gamma*(x-1) + delta*log(x)), col="red", add=T)  
alpha <- 1.03^(1:TT) * runif(TT,.95,1.05)                  
alpha <- 250 * alpha / alpha[1] / beta[1]

i <- rep(1:TT,TT:1); j <- sequence(TT:1); 
fi <- as.factor(i); fj <- as.factor(j)
mu.ij <- alpha[i] * beta[j]
Xij <- rpois(length(mu.ij), mu.ij)
xtabs(round(mu.ij)~i+j)        ## theoretical means
xtabs(Xij~i+j)                 ## generated run-off triangle 

## Chain Ladder (CL model)
CL <- glm(Xij~0+fi+fj, poisson)
exp(coef(CL))

# Question 7
beta.CL <- c(1, exp(coef(CL))[(TT+1):length(coef(CL))])
Hoerl <- glm(Xij~0+fi+I(j-1)+log(j), poisson)
coef(Hoerl)

# estimated parameters
gamma.Hoerl <- -coef(Hoerl)[TT+1] # 0.724162 
delta.Hoerl <- coef(Hoerl)[TT+2]  # 1.511653

# Question 8 

# (a) compute beta.Hoerl
k <- 1:TT
beta.Hoerl <- exp(-gamma.Hoerl*(k-1) + delta.Hoerl*log(k))

## (b) plot 
plot(beta.CL,  ylab = "Beta", main = "Hoerl curve")
curve(exp(-gamma*(x-1) + delta*log(x)), col="red", add=TRUE)
points(beta.Hoerl, pch = 19)
curve(exp(-gamma.Hoerl*(x-1) + delta.Hoerl*log(x)), 
      col="royalblue", add=TRUE)
legend("topright", c("Chain Ladder", "Hoerl"), 
       pch = c(1,19), cex = .75)

# Question 9
anova(Hoerl,CL)

# Question 10  
set.seed(1)

for(a in 1:10) {
    Xij <- rpois(length(mu.ij), mu.ij)
    CL <- glm(Xij~0+fi+fj, poisson)
    Hoerl <- glm(Xij~0+fi+I(j-1)+log(j), poisson)
    print(Hoerl$deviance - CL$deviance)
}

# Question 11
# (a)
Xij <- c(232,106,35,16,2,258,115,56,27,221,82,4,359,71,349)
i <- c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
j <- c(1,2,3,4,5,1,2,3,4,1,2,3,1,2,1)
TT <- 5

## create cumulative Xij
Xij <- as.numeric(unlist(tapply(Xij, i, cumsum)))

## generate the run-off triangle
fi <- as.factor(i); fj <- as.factor(j)
xtabs(Xij~i+j)

## Model 1 - Hoerl
Hoerl <- glm(Xij~0+fi+I(j-1)+log(j), poisson)

## Model 2 - Chain Ladder
CL <- glm(Xij~0+fi+fj, poisson)

## model comparison
anova(Hoerl,CL)
qchisq(.95,2)

## (b)
## model 11b - alpha.Hoerl
alpha.Hoerl <- glm(Xij~0+fj+I(i-1)+log(i), poisson)

## model comparison
anova(alpha.Hoerl,CL) 
qchisq(.95,2)

## Question 12
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
i <- rep(1:TT, TT:1); j <- sequence(TT:1); k <- i+j-1
fi <- as.factor(i); fj <- as.factor(j); fk <- as.factor(k)

## (a)
CL <- glm(Xij~fi+fj, poisson)
coefs.CL <- exp(coef(CL))
alpha.CL <- coefs.CL[1]*c(1, coefs.CL[2:TT])
beta.CL <- c(1, coefs.CL[(TT+1):length(coefs.CL)])


## (b)
round(alpha.CL %o% beta.CL)
## last column all zero                                

## (a) arithmetic separation
AS <- glm(Xij~fj+fk, poisson)

## (b) intercept term vs. observation
exp(coef(AS)[1])
Xij[1]                                          

cc <- exp(coef(AS))
beta.AS <- c(1,cc[2:8])*cc[1]; gamma.AS <- c(1,cc[9:15])
par(mfrow=c(1,2)); plot(gamma.AS); plot(log(gamma.AS))

ab <- coef(lm(log(gamma.AS)~I(1:8)))
par(mfrow=c(1,1)); plot(log(gamma.AS))
abline(ab[1],ab[2],col="red")

## Question 14
## (a)
l <- 1:15
gamma.fitted <- exp(ab[1]+l*ab[2])

## (b)
gammas <- c(gamma.AS, gamma.fitted[9:15])
jjj <- rep(1:8,8); kkk <- jjj + rep(1:8,each=8) - 1
mm <- beta.AS[jjj] * gammas[kkk] ## fitted values in vector form
round(matrix(mm,8,byrow=TRUE),3) ## fitted values in matrix form

## Question 15
## a)
Threeway <- glm(Xij~fi+fj+fk, poisson)

## b)
anova(CL,Threeway) # 12.907
qchisq(.95,6)      # [1] 12.59159

## c)
anova(AS,Threeway) # 18.958

## d)
AIC(CL, AS, Threeway)