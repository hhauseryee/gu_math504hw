setwd("C:/Users/yeeha/OneDrive/Documents/Graduate school/Math 504 - Numerical Methods/data")

## Problem 2

### 2 (c) Write an R function GramSchmidt(A) which returns the 
# matrix Q in (a). Check that your function works and com-
#  pare to the result of using R's qr function for some non-
# trivial choice of A.

GramSchmidt <- function(A){
  row_nums <- dim(A)[1]
  col_nums <- dim(A)[2]
  if (col_nums > row_nums){print("invalid matrix")}
  else{
    #browser()
    Q <- matrix(0, nrow = row_nums, ncol = col_nums)
    # populate first col of Q, j=1
    q_tilde <- A[,1]
    q <- q_tilde/sqrt(sum(q_tilde^2))
    Q[,1] <- q
    
    # j=2 and up
    for (j in 2:col_nums){
      q_tilde <- A[,j]
      for (i in 1:(j-1)){
        overlap <- - sum(A[,j]*Q[,i])*Q[,i]
        q_tilde <- q_tilde + overlap
      }
      q <- q_tilde/sqrt(sum(q_tilde^2))
      Q[,j] <- q
    }
    return(Q)
  }
}

A <- matrix(c(4,6,2,
              3,1,0,
              -1,3,-3,
              5,1,1),
            nrow = 4, ncol = 3, byrow = TRUE)

norm <- GramSchmidt(A)
norm

# check against R methods

qr_out <- qr(A)
Q_inR <- qr.Q(qr_out)
Q_inR


## Problem 4


# Attached are two files, senators_formatted.txt, which pro-
#  vides the names of all Senators in the 109th Senate of the U.S.
#along with their state and party affiliations (D=democrate,
# R=republican) and, votes_formatted.txt, which provides all
# votes for each senator over all bills considered (542). The
# 2nd through 101st column of votes_formatted.txt give the
# votes for each senator in the same order of senators as given in
# senators_formatted.txt, an entry of 1, 0, -1 corresponds to
# a YES, ABSENT, and NO votes. The first column gives the
# name of the bill voted on.

names <- read.csv("senators_formatted.txt", header = TRUE, sep = "",stringsAsFactors = FALSE)
votes <- read.csv("votes_formatted.txt",header = TRUE, sep = "")
senators <- t(votes)
# head(senators)

colnames(senators) <- senators[1,]
senators <- senators[-1, ] 

# assure that all entries are integers, not characters
mode(senators) <- "integer"

### 4. (a)

# center the data
mu <- rep(0,542)
for (i in 1:100){mu = mu + senators[i,]}
mu <- mu/100

# do a for loop to subtract mu one-by-one
for (i in 1:100){senators[i,] = senators[i,] - mu}

# calculate omega
omega <- matrix(0,542,542)
for (i in 1:100){omega = omega + senators[i,] %*% t(senators[i,])}

# calculate first two dominant eigenvectors
eigen(omega)$vectors[,1]
eigen(omega)$vectors[,2]

# corresponding evectors
evals <- eigen(omega)$values
evals[1:2]

# [1] 14974.90  2542.22

### 4. (b) Perform a 1-d PCA on the centered senator data. To do
# this, project each x(i) onto the dominant eigenvector and
# produce a 1-d plot of the senators. Color the senators
# according to party affiliation. What fraction of the total
# variance is captured by the PCA?

# pull the dominant eigenvector
evec1 <- eigen(omega)$vectors[,1]

# compute values for 1D plot
xp1 <- rep(0,100)
for (i in 1:100){xp1[i] <- sum(senators[i,] * evec1)}

# prepare data and colors
xp1 <- data.frame(xp1,1) ## 1 is your "height"

colors <- rep("c",100)
colors[names[,3] == "R"] <- "red"
colors[names[,3] == "D"] <- "blue"
colors[names[,3] == "I"] <- "green"

# 1D plot
plot(xp1, type="b", col=colors)

# What fraction of the total variance is captured by the PCA?

evals[1]/sum(evals)

### 4.   (c) Now repeat for a 2-d PCA and produce a 2-d plot.

# second most dominant eigenvector
evec2 <- eigen(omega)$vectors[,2]

# compute 2nd dimension values
xp2 <- rep(0,100)
for (i in 1:100){xp2[i] <- sum(senators[i,] * evec2)}

# the plot
plot(xp1[,1],xp2, col=colors)

# What fraction of the total variance is captured by the PCA?

sum(evals[1:2])/sum(evals)