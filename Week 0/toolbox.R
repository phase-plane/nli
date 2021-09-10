# Week 0 - A1

# read commands from a file
source(file = "test.R")

# to put the numbers written to a file (.txt OR csv) into a vector, use
vec <- scan(".../Desktop/anyfile.txt")

# identity matrix
diag(x=5, nrow=5, ncol=3)

# some matrix fun
vec <- 1:20
mat <- matrix(vec, ncol=5, byrow=TRUE); mat
ones <- matrix(1, nrow=5, ncol=4); ones

mat1 <- matrix(vec, ncol=5, byrow=TRUE) 
mat2 <- matrix(vec, ncol=5); mat2

diag(4) # 4x4 ID
diag(1:4) # dope
diag(nrow=3) # default is square matrix
diag(ones) <- 2; ones # superimpose diag. values on existing matrix

a <- 1:5
b <- 6:10
rbind(a,b)
cbind(a,b)

words <- c("Testing", "testing", "one", "two", "three")
for (p in words) {
  print(p)
}
?"for"

# fast and vectorised
x <- ifelse(cond, 1, 2)

# matrix operations
a.matrix <- rbind(10:15, 1:6)
a.matrix

# multiply matrices
b.matrix <- 2:7
a.matrix %*% b.matrix
# * element-wise
# a.matrix^2 (same)

# for illustration

# element-wise
1:3 * 1:3
t(1:3) * 1:3 
(1:3)^2      # [1] 1 4 9

# t() transpose

# inner product
1:3 %*% 1:3 
t(1:3) %*% 1:3
crossprod(1:3)

# but
crossprod(t(1:3)) 
1:3 %*% t(1:3)
1:3 %o% 1:3
tcrossprod(1:3)  # 3x3 matrix

words[c(3:5, 1)]

# select indices
numbers <- c(0, 3:5, 20, 0)
numbers[numbers<10]

a.matrix
a.matrix[1, numbers>4]

# testing equality 
all.equal(x,y)
isTRUE(all.equal())

# some arithmetic
signif(x, digits = 10)

# integration
integrate(dnorm, 0, +Inf)

set.seed(17); 
N <- rpois(10, lambda=4)  ## draw 10 Poisson(4) random numbers

# function computing the log-likelihood of mu with sample N
logLik <- function (mu) sum(dpois(N,mu,log=TRUE))

# find the maximum of this function
result <- optimize(f=logLik, lower=0, upper=1000, maximum=TRUE)
str(result)

test.fun <- function(x) sin(x)^2 + x/10
new.result <- optimise(f=test.fun, lower = 0, upper = 2*pi, maximum = FALSE)
str(new.result)
curve(sin(x)^2 + x/10, from = 0,to = 2*pi)

Y <- rnorm(1000)
hist(Y, breaks=seq(-4.0, 4.0, 0.5))