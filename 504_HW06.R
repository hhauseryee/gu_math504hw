setwd("C:/Users/yeeha/OneDrive/Documents/Graduate school/Math 504 - Numerical Methods")
library(knitr)
library(igraph)
library(Matrix)
library(readr)


## Problem 2

A <- matrix(c(1,2,5,2,9,7,5,7,4), nrow = 3, byrow = TRUE)

### 2. (a)
power_iter <- function(x0,A, tol=10^-8) {
  iter <- 0
  r1 <- 1 # initialize values for rayleigh quotients to some miscellaneous number
  r2 <- 2
  while (abs(r1-r2) > tol){
    r1 <- (t(x0) %*% A %*% x0)/sum(x0*x0) # Rayleigh quotient for x0
    x1 <- A %*% x0
    r2 <- (t(x1) %*% A %*% x1)/sum(x1*x1) # Rayleigh quotient for x1
    x0 <- x1/sqrt(sum(x1^2))
    iter <- iter +1
  }
  list(x0,iter)
}

de <- power_iter(c(1,1,0),A)
de

# This value of lambda was found very quickly, in 7 iterations.

# [,1]
# [1,] 0.3089619
# [2,] 0.7441196
# [3,] 0.5923078

# [[2]]
#[1] 9

### 2. (b)

evec <- unlist(de[1])
class(evec)

# compute the eigenvalue corresponding to the dominant eigenvector, using the Rayleigh quotient
evalue <- (t(evec) %*% A %*% evec)/(t(evec)%*%evec)
evalue

# [1,] 15.40229

### 2. (c)

eigen(A)
# $values
# [1] 15.402292  2.108321 -3.510612

# $vectors
# [,1]       [,2]       [,3]
# [1,] -0.3089629  0.6938644  0.6504569
# [2,] -0.7441200 -0.6022863  0.2890269
# [3,] -0.5923068  0.3947194 -0.7024025

# The approximation of the eigenvector using the power iteration and consqeuent computation of the e-value
# using the Rayleigh quotient are very close, identical up to the sixth decimal place.

## Problem 3

### 3. (a)


A <- matrix(c(1,1,0,1,1,
              0,1,1,0,0,
              1,0,1,1,0,
              0,1,1,1,1,
              1,0,0,0,1),
            byrow=TRUE, nrow = 5)
a <- base::rowSums(A)
A <- A/a
B <- matrix(rep(1/nrow(A),nrow(A)^2),
            byrow=TRUE, nrow = nrow(A))



# function to try several p-values
evector <- function(p,v,A,tol=10^-6){
  M <- (1-p)*A + p*(1/nrow(A))
  iter <- 0
  r1 <- 1
  r2 <- 2
  while (abs(r1-r2) > tol){
    r1 <- (t(v) %*% t(M) %*% v)/sum(v*v) # Rayleigh quotient for x0
    v <- t(M) %*% v
    r2 <- (t(v) %*% t(M) %*% v)/sum(v*v) # Rayleigh quotient for x1
    v <- v/sum(v)
    iter <- iter + 1
  }
  return(list(v, iter))
}

# test out some values of p
evector(.15,c(1,2,0,2,-1),A)

# [1,] 0.2205560
# [2,] 0.1977600
# [3,] 0.2105675
# [4,] 0.1733827
# [5,] 0.1977339

# [[2]]
# [1] 12

p <- .15
M <- (1-p)*A + p*B
eigen(t(M))

# [1] 1.0000000+0.0000000i 0.4250000+0.0000000i

evector(.2,c(1,2,0,2,-1),A)

# [,1]
# [1,] 0.2190269
# [2,] 0.1979592
# [3,] 0.2102355
# [4,] 0.1748506
# [5,] 0.1979278

# [[2]]
# [1] 11

p <- .2
M <- (1-p)*A + p*B
eigen(t(M))

# [1] 1.0000000+0.0000000i 0.4000000+0.0000000i 

evector(.25,c(1,2,0,2,-1),A)

# [1,] 0.2175294
# [2,] 0.1981587
# [3,] 0.2098631
# [4,] 0.1763312
# [5,] 0.1981175

# [1] 10

p <- .25
M <- (1-p)*A + p*B
eigen(t(M))

# [1] 1.0000000+0.0000000i 0.3750000+0.0000000i


# as p increases, the number of iterations needed to converge decreases, as the distance between them becomes greater and greater



# creating and testing a function to normalize the rows of a matrix
# sumnorm <- function(a){a/sum(a)}
# y <- matrix(c(1,5,7,2), byrow = TRUE, nrow =2)
# y <- apply(y, 1, sumnorm)


### 3. (b)

# read in data
hollins <- read_table2("data/edges.txt", col_names = FALSE)

# create the adjacency matrix
graph <- graph_from_data_frame(hollins) # what = c("edges", "vertices", "both")?
A1 <- as_adjacency_matrix(graph)

# add edges to A1, ensuring that no row is comprised of all zeros
for (i in 1:dim(A1)[1]){
  A1[i,i] <- 1
}

# normalize rows of adjacency matrix
v <- Matrix::rowSums(A1)
A1 <- A1/v


# power iteration function using M matrix
power_m <- function(p,v=rep(1,nrow(A)),A,tol=10^-6) {
  M <- (1-p)*A + p*(1/nrow(A))
  tM <- t(M)
  # iter <- 0
  for (i in 1:150){
    v <- tM %*% v
    v <- v/sum(v)
    # iter <- iter + 1
  }
  return(v)
}

# compute the dominant evector of t(M)
hollins_evec <- power_m(p=.15,A=A1)

which(hollins_evec[,1] == max(hollins_evec[,1]), arr.ind = TRUE)

largest <- sort(hollins_evec, decreasing = TRUE)
largest <- largest[1:5]
largest

for (i in largest) {
  maxes <- which(hollins_evec[,1] == i, arr.ind = TRUE)
  print(maxes)
}


# try out eigen()
tM <- t((1-.15)*A1 + .15*(1/nrow(A1)))
# eigen(tM)
# this causes my R to crash everytime it is run


## Problem 4

nd <- read.table("data/NotreDame.txt")
graph2 <- graph_from_data_frame(nd)

### 4.(a)

A2 <- as_adjacency_matrix(graph2, sparse = TRUE)


# add edges to A2, ensuring that no row is comprised of all zeros
length <- dim(A2)[1]
A2 <- cBind(A2,rep(1,length))
A2 <- rBind(A2,rep(1,length+1))
dim(A2)
class(A2)

# normalize rows of adjacency matrix
v2 <- Matrix::rowSums(A2)
A2 <- A2/v2

# A2 <- as_adjacency_matrix(graph2, sparse = FALSE)
# yields a fatal error which causes RStudio to crash

### 4. (b)

# for large, computationally intensive M

power_M <- function(p,v=rep(1,nrow(A)),A,tol=10^-6) {
  # iter <- 0
  # browser()
  for (i in 1:100){
    bv <- sum(p*(1/nrow(A))*v)
    v <- (1-p)*t(A) %*% v + rep(bv,nrow(A))
    v <- v/sum(v)
    # iter <- iter + 1
  }
  return(v)
}

# v <- (1-p)*A %*% v + p*(1/nrow(A))*v

nd_vec <- power_M(p=.15,A=A2)
nd_vec

# find the index of the most popular pages
which(nd_vec[,1] == max(nd_vec[,1]), arr.ind = TRUE)

# 325730

# Since it isn't giving node IDs, tried this to correct the situation, but names() doesn't work for sparse matrices
# names(nd_vec) <- colnames(A2)

# finds indices of 4 most popular pages
top_four <- nd_vec[nd_vec[,1] > .0019]
top_four

for (i in top_four){
  index <- which(nd_vec[,1] == i, arr.ind = TRUE)
  print(index)
}


