# generating a pseudo-random insurance portfolio

n.obs <- 10000; set.seed(4)

# class for gender
sx <- sample(1:2, n.obs, repl=TRUE, prob=c(6,4))
sx <- as.factor(sx)

# class for job
jb <- as.factor(sample(1:3, n.obs, replace = TRUE, prob = c(3,2,1)))

# factors with correlation (region, type of car)
re.tp <- sample(1:9, n.obs, repl=TRUE, 
                prob=c(.1,.05,.15,.15,.1,.05,.1,.1,.2))

tp <- c(1,2,3,1,2,3,1,2,3)[re.tp]; tp <- as.factor(tp)                # *
re <- c(1,1,1,2,2,2,3,3,3)[re.tp]; re <- as.factor(re)

# view our generated content (order matters (x,y))
table(list(region=re, type=tp))

# check memory usage
object.size(re.tp)
rm(re.tp)

# monthly duration of policy
mo <- 3 * sample(1:4, n.obs, repl=TRUE, prob=c(1,1,0,8))

# stated facts to construct mu for our simulation                     # **
mu <- 0.05 * c(1,1.2)[sx] *
  c(1,1,1)[jb] *
  c(1,1.2,1.44)[re] *
  1.2^(0:2)[tp] *
  mo/12

# table of the number of claims per policy
# the number of claims in each cell is a poisson process
y <- rpois(n.obs, mu) 
table(y)
var(y)/mean(y)
table(list(numCl=y,gender=sx))

# interaction operator ':'
table(list(numCl=y,gender.region=sx:re))
table(list(numCl=y,gender.region=interaction(sx,re)))
