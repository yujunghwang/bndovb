#' @title bndovb
#' @description This function runs a two sample least squares when auxiliary data contains every right-hand side regressor
#' and main data contains a dependent variable and every right-hand side regressor but one omitted variable.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references Hwang, Yujung (2021). Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.
#' @importFrom utils install.packages
#' @import stats
#' @import np
#' @import pracma
#' @import dplyr
#' @import MASS
#'
#' @param maindat Main data set
#' @param auxdat Auxiliary data set
#' @param depvar A name of a dependent variable in main dataset
#' @param ovar A name of an omitted variable in main dataset which exists in auxiliary data
#' @param comvar A vector of common regressors existing in both main data and auxiliary data
#' @param method CDF and Quantile function estimation method.
#' Users can choose either 1 or 2. If the method is 1, the CDF and quantile function is estimated assuming a parametric normal distribution.
#' If the method is 2, the CDF and quantile function is estimated using a nonparaemtric estimator in Li and Racine(2008), Li, Lin, and Racine(2013).
#' Default is 1.
#'
#' @return Returns a list of 2 components : \describe{
#' \item{hat_beta_l}{lower bound estimates of regression coefficients}
#'
#' \item{hat_beta_u}{upper bound estimates of regression coefficients}}
#'
#' @export
bndovb <- function(maindat,auxdat,depvar,ovar,comvar,method=1){

  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("np")
  requireNamespace("pracma")

  #############
  # check if inputs are there in a correct form
  #############

  if (!is.matrix(maindat) & !is.data.frame(maindat)){
    stop("please provide main data in either matrix or data frame format.")
  }

  if (!is.matrix(auxdat) & !is.data.frame(auxdat)){
    stop("please provide auxiliary data in either matrix or data frame format.")
  }

  # check if column names of auxiliary data exists
  if (is.null(colnames(auxdat))){
    stop("column names of auxiliary data do not exist.")
  }

  # check if column names of main data exists
  if (is.null(colnames(maindat))){
    stop("column names of main data do not exist.")
  }

  # check if auxiliary dataset includes every independent regressor
  if ((sum(comvar%in%colnames(auxdat))<length(comvar)) | !(ovar%in%colnames(auxdat)) ){
    stop("auxiliary dataset does not contain every right-hand side regressor.")
  }

  # check if main dataset includes every independent regressor
  if (sum(comvar%in%colnames(maindat))<length(comvar)){
    stop("main dataset does not contain every common right-hand side regressor.")
  }

  # check if main dataset includes dependent variable
  if (!(depvar%in%colnames(maindat))){
    stop("main dataset does not include the dependent variable.")
  }

  # check if method is specified correctly
  if (!(method%in%c(1,2))){
    stop("Incorrect method was specified. Method should be either 1 or 2.")
  }

  #############
  # prepare data in a right form
  #############

  # number of observations
  Nm <- dim(maindat)[1]
  Na <- dim(auxdat)[1]

  # add 1 vector
  comvar <- c(comvar,"con")
  maindat$con <- rep(1,Nm)
  auxdat$con <- rep(1,Na)

  # leave only necessary variables and make the order of variables consistent
  maindat <- maindat[,c(depvar,comvar)]
  auxdat <- auxdat[,c(ovar,comvar)]

  # number of regressors in a regrssion model
  nr <- length(comvar)+length(ovar)

  #############
  # estimate CDF and Quantile function
  #############

  if (method==1){

    # estimate N(depvar | comvar)
    f1 <- paste0(depvar,"~ 0 +",comvar[1])
    if (length(comvar)>1){
      for (k in 2:length(comvar)){
        f1 <- paste0(f1,"+",comvar[k])
      }
    }
    oout1 <- lm(formula=f1,data=maindat) ## regression without intercept because of "con" in "comvar"
    Fypar <- matrix(oout1$coefficients,ncol=1)
    yhat  <- as.matrix(maindat[,comvar])%*%Fypar
    ysd   <- sd(oout1$residuals,na.rm=TRUE)

    # estimate N(ovar | comvar)
    f2 <- paste0(ovar,"~ 0 +",comvar[1])
    if (length(comvar)>1){
        for (k in 2:length(comvar)){
          f2 <- paste0(f2,"+",comvar[k])
        }
    }
    oout2 <- lm(formula=f2,data=auxdat) ## regression without intercept because of "con" in "comvar"
    Fopar <- matrix(oout2$coefficients,ncol=1)

    # prediction in main data, not auxiliary data
    ohat  <- as.matrix(maindat[,comvar])%*%Fopar
    osd   <- sd(oout2$residuals,na.rm=TRUE)

    #############
    # compute bounds of E[(depvar)*(omitted variable)]
    #############

    ovar_m_l <- rep(NA,Nm)
    ovar_m_u <- rep(NA,Nm)

    for (k in 1:Nm){
      ovar_m_u[k] <- qnorm(p=   pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd) ,mean=ohat[k],sd=osd)
      ovar_m_l[k] <- qnorm(p=(1-pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd)),mean=ohat[k],sd=osd)
    }

  } else if (method==2){

    # estimate f(depvar | comvar) nonparametrically
    # bandwidth selection
    bws1 <- npcdistbw(ydat=maindat[,depvar],xdat=maindat[,comvar])
    Fyz  <- npcdist(bws1)$condist ### Fyz$condist saves the predicted cdf values

    bws2 <- npcdistbw(ydat=auxdat[,ovar],xdat=auxdat[,comvar])

    # compute matching function mu(depvar) = ovar| depvar, comvar
    mu_y <- function(xx,ccdf,maximize){
      if (maximize==1){
        # find matching ovar to maximize E[depvar * ovar]
        ovar1  <- npqreg(bws2,tau=ccdf,exdat=xx)$quantile
      } else {
        # find matching ovar to minimize E[depvar * ovar]
        ovar1  <- npqreg(bws2,tau=(1-ccdf),exdat=xx)$quantile
      }
      return(ovar1)
    }

    ovar_m_l <- rep(NA,Nm)
    ovar_m_u <- rep(NA,Nm)

    for(i in 1:Nm){
      eexdat <- data.frame(maindat[i,comvar])
      colnames(eexdat) <- c(1:length(comvar)) ### make the data frame similar to txdat

      ovar_m_l[i]     <- mu_y(eexdat,ccdf=Fyz[i],maximize=0)
      ovar_m_u[i]     <- mu_y(eexdat,ccdf=Fyz[i],maximize=1)
      rm(eexdat)
    }

  } else {
    stop("Method should be either 1 or 2.")
  }

  #############
  # compute lower bound and upper bound
  #############

  mu_l <- mean(maindat[,depvar]*ovar_m_l,na.rm=TRUE)
  mu_u <- mean(maindat[,depvar]*ovar_m_u,na.rm=TRUE)

  hat_beta_l <- rep(NA,nr)
  hat_beta_u <- rep(NA,nr)

  # submatrices
  A1 <- (t(as.matrix(auxdat[,ovar]))%*%as.matrix(auxdat[,ovar]))/Na
  A2 <- (t(as.matrix(auxdat[,ovar]))%*%as.matrix(auxdat[,comvar]))/Na
  C  <- as.matrix(rbind(maindat[,comvar],auxdat[,comvar]))
  A3 <- (t(C)%*%C)/(Na+Nm)

  XX <- as.matrix(rbind(cbind(A1,A2),cbind(t(A2),A3)))

  # OLS formula
  B <- (t(as.matrix(maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/Nm
  B_l <- matrix(c(mu_l,B),ncol=1)
  B_u <- matrix(c(mu_u,B),ncol=1)

  hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
  hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)

  colnames(hat_beta_l) <- c(ovar,comvar)
  colnames(hat_beta_u) <- c(ovar,comvar)

  # change the order of OLS coefficients
  comvar2 <- comvar[comvar!="con"]
  hat_beta_l <- c(hat_beta_l[,"con"],hat_beta_l[,ovar],hat_beta_l[,comvar2])
  hat_beta_u <- c(hat_beta_u[,"con"],hat_beta_u[,ovar],hat_beta_u[,comvar2])

  return(list(hat_beta_l=hat_beta_l,hat_beta_u=hat_beta_u))
}
