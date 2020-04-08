setwd("C:/Users/yeeha/OneDrive/Documents/Graduate school/Math 504 - Numerical Methods")

# Import and set up training data
mtrain <- as.matrix(read.csv("mnist_train.csv", header=F))
number <- mtrain[,1]
X <- mtrain[,2:ncol(mtrain)]
colnames(X) <- NULL

# add some noise to "solve" the invertibility problem
set.seed(585)
X <- X + matrix(.1*rnorm(nrow(X)*ncol(X)), nrow=nrow(X))
ones <- rep(1,60000)
X <- cbind(ones,X)
dim(X)

# set up y
number <- mtrain[,1]
y <- 1*(number == 3)

# set up xs
x <- X[13,]


#log L(a)
logL <- function(x,y,a){
  llikely <- 0
  for (i in 1:60000){
    arg <- -sum(a*x[i,])
    llikely <- llikely + (1-y[i])*arg - log (1 +exp(arg))
  }
  return(llikely)
}

# logL(X,y,a=rep(0,785))
# system.time(logL(X,y,a=rep(0,785)))

#user  system elapsed 
#3.29    0.03    1.06 

#Gradient of log L(a)
grad_logL <- function(x,y,a=rep(0,785)) {
  g <- matrix(0, nrow=nrow(x), ncol = ncol(x))
  for (i in 1:60000){
    arg <- exp(-sum(a*x[i,]))
    g[i,] <- -x[i,]*(1-y[i]) + x[i,]*(arg/(1+arg))
  }
  return(colSums(g))
}

# grad_logL(X,y)
# system.time(grad_logL(X,y))
# user  system elapsed 
# 3.55    0.03    3.83


hessian_logL <- function(x,a=rep(0,785)){
  H <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  for (i in 1:60000){
    arg <- exp(sum(-a*x[i,]))
    H <- H - (x[i,] %*% t(x[i,])) *(arg/(1+arg)^2)
  }
  return(H)
}

# H <- hessian_logL(X)  
# system.time(hessian_logL(X))

# user  system elapsed 
# 158.42    8.76  177.44 

# damped Newton's method for optimization
damped_N <- function(x,y,a0=rep(0,785)) {
  # browser()
  for (i in 1:11){
    g <- grad_logL(x,y,a0)
    H <- hessian_logL(x,a0)
    d <- - solve(H,g)
    s <- 1
    while (logL(X,y,a0 + s*d) < logL(X,y,a0)){
      s <- s/2
    }
    a1 <- a0 + s*d
    a0 <- a1
    }
  return(a1)
}

# with 12 iterations
# opt_a <- damped_N(X,y,rep(0,785))

# with 11 iterations
opt_11 <- damped_N(X,y,rep(0,785))

logL(X,y,opt_11)

# -4357.228

logL(X,y,rep(0,785))

# -41588.83

### 4. (d)

classifier <- function(x,p,a){
  classified <- rep(0,nrow(X))
  for (i in 1:nrow(X)) {
    classified[i] <- 1* (1/(1 + exp(sum(-a*x[i,]))) >= p)
  }
  return(classified)
}

test <- classifier(X,.3,opt_11)

accuracy <- function(classifications,y){
  correct <- sum(1*(classifications == y))
  return(correct/length(y))
}

accuracy(test,y)
