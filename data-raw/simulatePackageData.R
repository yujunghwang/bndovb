# simulate package data

library(MASS)
library(usethis)

set.seed(210601)

# sample size
Nm <- 100000 # main data
Na <- 50000 # auxiliary data

## Example data for a function 'bndovb'
# use same DGP in maindat and auxdat
maindat_nome <- mvrnorm(Nm,mu=c(2,3,1),Sigma=rbind(c(2,1,1),c(1,2,1),c(1,1,2)))
auxdat_nome  <- mvrnorm(Na,mu=c(2,3,1),Sigma=rbind(c(2,1,1),c(1,2,1),c(1,1,2)))

maindat_nome <- as.data.frame(maindat_nome)
auxdat_nome <- as.data.frame(auxdat_nome)

colnames(maindat_nome) <- c("x1","x2","x3")
colnames(auxdat_nome) <- c("x1","x2","x3")

# this is a true parameter which we try to get bounds on
truebeta <- matrix(c(2,1,3,2),ncol=1)

# generate a dependent variable
maindat_nome$y <- as.matrix(cbind(rep(1,Nm),maindat_nome[,c("x1","x2","x3")]))%*%truebeta

# main data misses one omitted variable "x1"
maindat_nome <- maindat_nome[,c("x2","x3","y")]

use_data(maindat_nome,auxdat_nome,overwrite=TRUE)


## Example data for a function 'bndovbme', continuous proxy variables
# set DGP
nu      <- 0.5   # sd of measurement errors in proxy variables
beta    <- c(0,1,1,1) # true parameters in a regression model
gamma   <- c(0,1,1) # parameters to generate correlation between covariates
samsize <- c(6000)  # sample size
mu      <- c(0,0,0,0) # average of covariates
sigma   <- eye(4)

#### simulate data
A <- rbind( c(1,0,0,0), c(0,1,0,0), c(gamma[2],gamma[3],1,0), c(beta[3]+beta[2]*gamma[2],beta[4]+beta[2]*gamma[3], beta[2],1))
B <- c(0,0,gamma[1],beta[1])
mu2    <- A%*%mu + B
sigma2 <- A%*%sigma%*%t(A)
Sim     = 100         ;  # number of Monte Carlo simulations

n=6000;na=3000;nb=3000
simdata <- mvrnorm(n,mu=mu2,Sigma=sigma2)

w1<-simdata[,1]
w2<-simdata[,2]
x <-simdata[,3]
y <-simdata[,4]

# main data
w1_a <- w1[1:na]
w2_a <- w2[1:na]
x_a  <- x[ 1:na]
y_a  <- y[ 1:na]

# auxiliary data
w1_b <- w1[(na+1):n]
w2_b <- w2[(na+1):n]
x_b  <- x[ (na+1):n]
y_b  <- y[ (na+1):n]

# generate continuous proxies
z_b <- w2_b + cbind(rnorm(n-na,mean=0,sd=nu), rnorm(n-na,mean=0,sd=nu), rnorm(n-na,mean=0,sd=nu))

# main data does not include a variable w2
maindat_mecont <- data.frame(y=y_a,x=x_a,w1=w1_a)

# auxiliary data does not include a dependent variable y
# auxiliary data contain three proxy variables for the omitted variable w2
auxdat_mecont <- data.frame(x=x_b,w1=w1_b,z1=z_b[,1],z2=z_b[,2],z3=z_b[,3])

# save as a package data
use_data(maindat_mecont,auxdat_mecont,overwrite=TRUE)


## Example data for a function 'bndovbme', discrete proxy variables
# set DGP

n=6000;na=3000;nb=3000
mu2 <- c(0,0,0)
Sigma2 <- rbind(c(1,0.5,0.5),c(0.5,1,0.5),c(0.5,0.5,1))

simdata <- mvrnorm(n,mu=mu2,Sigma=Sigma2)
# discretize
simdata[,2] <- (simdata[,2]>0)+1

beta    <- c(0,1,1,1) # true parameters to get bounds on

# simulate a dependent variable
y <- cbind(rep(1,n),simdata)%*%as.matrix(beta) + rnorm(n)

w1<-simdata[,1]
w2<-simdata[,2]
x <-simdata[,3]

# main data
w1_a <- w1[1:na]
w2_a <- w2[1:na]
x_a  <- x[ 1:na]
y_a  <- y[ 1:na]

# auxiliary data
w1_b <- w1[(na+1):n]
w2_b <- w2[(na+1):n]
x_b  <- x[ (na+1):n]
y_b  <- y[ (na+1):n]

# set measurement matrices for discrete proxy variables
M_param <- list()
M_param[[1]] <- rbind(c(0.9,0.1),c(0.1,0.9))
M_param[[2]] <- rbind(c(0.9,0.1),c(0.1,0.9))
M_param[[3]] <- rbind(c(0.9,0.1),c(0.1,0.9))

CM_param <- list()
CM_param[[1]] <- t(apply(M_param[[1]],1,cumsum))
CM_param[[2]] <- t(apply(M_param[[2]],1,cumsum))
CM_param[[3]] <- t(apply(M_param[[3]],1,cumsum))

# simulate proxy variables
z_b <- matrix(NA,nrow=nb,ncol=3)
for (k in 1:nb){
  for (l in 1:3){
    z_b[k,l] <- which(runif(1)<CM_param[[l]][w2_b[k],])[1]
  }
}

# main data
maindat_medisc <- data.frame(y=y_a,x=x_a,w1=w1_a)
# auxiliary data
auxdat_medisc <- data.frame(x=x_b,w1=w1_b,z1=z_b[,1],z2=z_b[,2],z3=z_b[,3])

# save as a package data
use_data(maindat_medisc,auxdat_medisc)
