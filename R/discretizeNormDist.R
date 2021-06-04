#' @title discretizeNormDist
#' @description This function discretizes a normal distribution into n grid points with equal probability
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Adda, J., and Cooper, R. W. (2003)}{Dynamic economics: quantitative methods and applications. MIT press.}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.}}
#' @importFrom utils install.packages
#' @import stats
#'
#' @param n Number of grid points
#' @param mean Mean of a normal distribution to discretize
#' @param var Variance of a normal distribution to discretize
#' 
#' @return Returns an equal probability grid point : \describe{
#' \item{EVbetweenZ}{n grid points with equal probability}}
#'
#' @examples
#' discretizeNormDist(n=3,mean=0,var=1)
#' 
#' @export
discretizeNormDist <- function(n,mean,var){
  # discretize normal distribution into n grid points with equal prob
  
  Z <- rep(NA,n+1)
  EVbetweenZ <- rep(NA,n)
  
  # cutoff points
  Z[1]   <- -Inf
  Z[n+1] <-  Inf
  
  # Now recursively get the rest of the points
  for (ixi in 2:n){
    Z[ixi] <- sqrt(var) * qnorm((ixi-1)/n) + mean
  }
  stdZ <- (Z-mean)/sqrt(var)
  
  # Finding the expected value within each interval (see Adda & Cooper page 58)
  PDF <- dnorm(stdZ)
  for (ixi in 1:n){
    EVbetweenZ[ixi] <- n*sqrt(var)*(PDF[ixi]-PDF[ixi+1]) + mean
  }
  
  return(EVbetweenZ)
}