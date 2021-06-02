#' @title bndovb
#' @description This function runs a two sample least squares when auxiliary data contains every right-hand side regressor
#' and main data contains a dependent variable and every right-hand side regressor but one omitted variable.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.}}
#' @importFrom utils install.packages
#' @import stats
#' @import np
#' @importFrom pracma pinv eye
#' @importFrom MASS mvrnorm
#'
#' @param maindat Main data set
#' @param auxdat Auxiliary data set
#' @param depvar A name of a dependent variable in main dataset
#' @param ovar A name of an omitted variable in main dataset which exists in auxiliary data
#' @param comvar A vector of common regressors existing in both main data and auxiliary data
#' @param method CDF and Quantile function estimation method.
#' Users can choose either 1 or 2. If the method is 1, the CDF and quantile function is estimated assuming a parametric normal distribution.
#' If the method is 2, the CDF and quantile function is estimated using a nonparaemtric estimator in Li and Racine(2008) \doi{10.1198/073500107000000250}, Li, Lin, and Racine(2013) \doi{10.1080/07350015.2012.738955}.
#' Default is 1.
#' @param mainweights An optional weight vector for the main dataset. The vector length must be equal to the number of rows of 'maindat'
#' @param auxweights An optional weight vector for the auxiliary dataset. The vector length must be equal to the number of rows of 'auxdat'
#'
#' @return Returns a list of 2 components : \describe{
#' \item{hat_beta_l}{lower bound estimates of regression coefficients}
#'
#' \item{hat_beta_u}{upper bound estimates of regression coefficients}}
#'
#' @examples
#' data(maindat_nome)
#' data(auxdat_nome)
#'
#' bndovb(maindat=maindat_nome,auxdat=auxdat_nome,depvar="y",ovar="x1",comvar=c("x2","x3"),method=1)
#'
#'
#' @export
bndovb <- function(maindat,auxdat,depvar,ovar,comvar,method=1,mainweights=NULL,auxweights=NULL){

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

  if (length(ovar)>1){
    stop("there are too many omitted variables.")
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

  if (!is.null(mainweights)){
    # check if the weight vector for main data has a correct length
    if (length(mainweights)!=dim(maindat)[1]){
     stop("Incorrect length for the main data weight vector. The length must be equal to the number of rows of 'maindat'.")
    }
    # check if any weight vector includes NA or NaN or Inf
    if (sum(is.na(mainweights))>0|sum(is.nan(mainweights))>0|sum(is.infinite(mainweights))>0){
      stop("mainweights vector can not include any NAs or NaNs or Infs.")
    }
  }

  if (!is.null(auxweights)){
    # check if the weight vector for auxiliary data has a correct length
    if (length(auxweights)!=dim(auxdat)[1]){
      stop("Incorrect length for the auxiliary data weight vector. The length must be equal to the number of rows of 'auxdat'.")
    }
    # check if any weight vector includes NA or NaN or Inf
    if (sum(is.na(auxweights))>0|sum(is.nan(auxweights))>0|sum(is.infinite(auxweights))>0){
      stop("auxweights vector can not include any NAs or NaNs or Infs.")
    }
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
    if (is.null(mainweights)){
      oout1 <- lm(formula=f1,data=maindat) ## regression without intercept because of "con" in "comvar"
    } else{
      oout1 <- lm(formula=f1,data=maindat,weights=mainweights) ## regression without intercept because of "con" in "comvar"
    }

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
    if (is.null(auxweights)){
      oout2 <- lm(formula=f2,data=auxdat) ## regression without intercept because of "con" in "comvar"
    } else{
      oout2 <- lm(formula=f2,data=auxdat,weights=auxweights) ## regression without intercept because of "con" in "comvar"
    }
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
      if (!is.na(maindat[k,depvar]) & !is.nan(maindat[k,depvar]) & !is.na(yhat[k]) & !is.nan(yhat[k]) & !is.na(ysd) & !is.nan(ysd) & !is.na(ohat[k]) & !is.nan(ohat[k]) & !is.na(osd) & !is.nan(osd) ){
        ovar_m_u[k] <- qnorm(p=   pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd) ,mean=ohat[k],sd=osd)
        ovar_m_l[k] <- qnorm(p=(1-pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd)),mean=ohat[k],sd=osd)
      }
    }

  } else if (method==2){

    ### use np package

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

  # replace missing values to 0 and create a dummy for missingness
  Imaindat <- !is.na(maindat)
  Iauxdat  <- !is.na(auxdat)

  colnames(Imaindat) <- colnames(maindat)
  colnames(Iauxdat)  <- colnames(auxdat)

  maindat[!Imaindat] <-0
  auxdat[!Iauxdat]   <-0

  Iovar_m_l <- !is.na(ovar_m_l)
  Iovar_m_u <- !is.na(ovar_m_u)

  ovar_m_l[!Iovar_m_l] <-0
  ovar_m_u[!Iovar_m_u] <-0

  hat_beta_l <- rep(NA,nr)
  hat_beta_u <- rep(NA,nr)

  if (is.null(mainweights)){

    mu_l <- sum(maindat[,depvar]*ovar_m_l) / sum(Imaindat[,depvar]*Iovar_m_l)
    mu_u <- sum(maindat[,depvar]*ovar_m_u) / sum(Imaindat[,depvar]*Iovar_m_u)

  } else{

    mu_l <- sum(maindat[,depvar]*ovar_m_l*mainweights) / sum(Imaindat[,depvar]*Iovar_m_l*mainweights)
    mu_u <- sum(maindat[,depvar]*ovar_m_u*mainweights) / sum(Imaindat[,depvar]*Iovar_m_u*mainweights)

  }

  # submatrices
  if (is.null(auxweights)){

    A1 <- (t(as.matrix(auxdat[,ovar]))%*%as.matrix(auxdat[,ovar]))  /(t(as.matrix(Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,ovar]))
    A2 <- (t(as.matrix(auxdat[,ovar]))%*%as.matrix(auxdat[,comvar]))/(t(as.matrix(Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,comvar]))

  } else{

    A1 <- (t(as.matrix(auxweights*auxdat[,ovar]))%*%as.matrix(auxdat[,ovar]))  /sum(t(as.matrix(auxweights*Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,ovar]))
    A2 <- (t(as.matrix(auxweights*auxdat[,ovar]))%*%as.matrix(auxdat[,comvar]))/sum(t(as.matrix(auxweights*Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,comvar]))

  }


  if (is.null(auxweights) & is.null(mainweights)){

    C  <- as.matrix(rbind( maindat[,comvar], auxdat[,comvar]))
    IC <- as.matrix(rbind(Imaindat[,comvar],Iauxdat[,comvar]))

    A3 <- (t(C)%*%C)/(t(IC)%*%IC)

  } else if(!is.null(auxweights) & is.null(mainweights)){

    aw <- matrix(rep(auxweights, length(comvar)),ncol=length(comvar)) *(1/sum(auxweights)) * Na

    C  <- as.matrix(rbind( maindat[,comvar],aw* auxdat[,comvar]))
    IC <- as.matrix(rbind(Imaindat[,comvar],aw*Iauxdat[,comvar]))

    A3 <- (t(C)%*%C)/(t(IC)%*%IC)

  } else if(is.null(auxweights) & !is.null(mainweights)){

    mw <- matrix(rep(mainweights,length(comvar)),ncol=length(comvar)) *(1/sum(mainweights)) * Nm

    C  <- as.matrix(rbind(mw* maindat[,comvar],  auxdat[,comvar]))
    IC <- as.matrix(rbind(mw*Imaindat[,comvar], Iauxdat[,comvar]))

    A3 <- (t(C)%*%C)/(t(IC)%*%IC)

  } else{

    mw <- matrix(rep(mainweights,length(comvar)),ncol=length(comvar)) *(1/sum(mainweights)) * Nm
    aw <- matrix(rep(auxweights, length(comvar)),ncol=length(comvar)) *(1/sum(auxweights))  * Na

    C  <- as.matrix(rbind(mw* maindat[,comvar], aw* auxdat[,comvar]))
    IC <- as.matrix(rbind(mw*Imaindat[,comvar], aw*Iauxdat[,comvar]))

    A3 <- (t(C)%*%C)/(t(IC)%*%IC)

  }

  XX <- as.matrix(rbind(cbind(A1,A2),cbind(t(A2),A3)))

  # OLS formula
  if (is.null(mainweights)){
    B <- (t(as.matrix(maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/(t(as.matrix(Imaindat[,depvar]))%*%as.matrix(Imaindat[,comvar]))
  } else{
    B <- (t(as.matrix(mainweights*maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/(t(as.matrix(mainweights*Imaindat[,depvar]))%*%as.matrix(Imaindat[,comvar]))
  }

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
