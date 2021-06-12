#' @title bndovbme
#' @description This function runs a two sample least squares when main data contains a dependent variable and
#' every right hand side regressor but one omitted variable.
#' The function requires an auxiliary data which includes every right hand side regressor but one omitted variable,
#' and enough proxy variables for the omitted variable.
#' When the omitted variable is continuous, the auxiliary data must contain at least two continuous proxy variables.
#' When the omitted variable is discrete, the auxiliary data must contain at least three continuous proxy variables.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.}}
#' @importFrom utils install.packages
#' @import stats
#' @importFrom pracma pinv eye
#' @importFrom MASS mvrnorm
#' @import factormodel
#' @importFrom nnet multinom
#'
#' @param maindat Main data set. It must be a data frame.
#' @param auxdat Auxiliary data set. It must be a data frame.
#' @param depvar A name of a dependent variable in main dataset
#' @param pvar A vector of the names of the proxy variables for the omitted variable.
#' When proxy variables are continuous, the first proxy variable is used as an anchoring variable.
#' When proxy variables are discrete, the first proxy variable is used for initialization (For details, see a documentation for "dproxyme" function).
#' @param ptype Either 1 (continuous) or 2 (discrete). Whether proxy variables are continuous or discrete. Default is 1 (continuous).
#' @param comvar A vector of the names of the common regressors existing in both main data and auxiliary data
#' @param sbar A cardinality of the support of the discrete proxy variables. Default is 2. If proxy variables are continuous, this variable is irrelevant.
#' @param mainweights An optional weight vector for the main dataset. The length must be equal to the number of rows of 'maindat'.
#' @param auxweights An optional weight vector for the auxiliary dataset. The length must be equal to the number of rows of 'auxdat'.
#' @param normalize Whether to normalize the omitted variable to have mean 0 and standard deviation 1. Set TRUE or FALSE.
#' Default is TRUE. If FALSE, then the scale of the omitted variable is anchored with the first proxy variable in pvar list.
#' @param signres An option to impose a sign restriction on a coefficient of an omitted variable. Set either NULL or pos or neg.
#' Default is NULL. If NULL, there is no sign restriction.
#' If 'pos', the estimator imposes an extra restriction that the coefficient of an omitted variable must be positive.
#' If 'neg', the estimator imposes an extra restriction that the coefficient of an omitted variable must be negative.
#'
#' @return Returns a list of 4 components : \describe{
#' \item{hat_beta_l}{lower bound estimates of regression coefficients}
#'
#' \item{hat_beta_u}{upper bound estimates of regression coefficients}
#'
#' \item{mu_l}{lower bound estimate of E\[ovar*depvar\]}
#'
#' \item{mu_u}{upper bound estimate of E\[ovar*depvar\]}}
#'
#' @examples
#' ## load example data
#' data(maindat_mecont)
#' data(auxdat_mecont)
#'
#' ## set ptype=1 for continuous proxy variables
#'  pvar<-c("z1","z2","z3")
#'  cvar<-c("x","w1")
#' bndovbme(maindat=maindat_mecont,auxdat=auxdat_mecont,depvar="y",pvar=pvar,ptype=1,comvar=cvar)
#'
#' ## set ptype=2 for discrete proxy variables
#' data(maindat_medisc)
#' data(auxdat_medisc)
#' bndovbme(maindat=maindat_medisc,auxdat=auxdat_medisc,depvar="y",pvar=pvar,ptype=2,comvar=cvar)
#'
#' @export
bndovbme <- function(maindat,auxdat,depvar,pvar,ptype=1,comvar,sbar=2,mainweights=NULL,auxweights=NULL,normalize=TRUE,signres=NULL){

  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("pracma")
  requireNamespace("factormodel")
  requireNamespace("nnet")

  #############
  # check if inputs are there in a correct form
  #############

  if (!is.data.frame(maindat)){
    stop("please provide main data in a data frame format.")
  }

  if (!is.data.frame(auxdat)){
    stop("please provide auxiliary data in a data frame format.")
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
  if ((sum(comvar%in%colnames(auxdat))<length(comvar)) | (sum(pvar%in%colnames(auxdat))<length(pvar)) ){
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

  # check if the proxy variable type is correctly specified
  if (!(ptype%in%c(1,2))){
    stop("Incorrect type was specified for proxy variables. ptype should be either 1 or 2.")
  }

  # check if there are enough proxy variables
  if ((ptype==1) & (length(pvar)<2)){
    stop("There are insufficient number of proxy variables. There must be at least 2 proxy variables when the omitted variable is continuous.")
  }

  if ((ptype==2) & (length(pvar)<3)){
    stop("There are insufficient number of proxy variables. There must be at least 3 proxy variables when the omitted variable is discrete.")
  }

  if (!is.null(mainweights)){
    # check if the weight vector has right length
    if (length(mainweights)!=dim(maindat)[1]){
      stop("The length of 'mainweights' is not equal to the number of rows of 'maindat'.")
    }
    # check if any weight vector includes NA or NaN or Inf
    if (sum(is.na(mainweights))>0|sum(is.nan(mainweights))>0|sum(is.infinite(mainweights))>0){
      stop("mainweights vector can not include any NAs or NaNs or Infs.")
    }
  }

  if (!is.null(auxweights)){
    # check if the weight variable is included in the auxdat
    if (length(auxweights)!=dim(auxdat)[1]){
      stop("The length of 'auxweights' is not equal to the number of rows of 'auxdat'.")
    }
    # check if any weight vector includes NA or NaN or Inf
    if (sum(is.na(auxweights))>0|sum(is.nan(auxweights))>0|sum(is.infinite(auxweights))>0){
      stop("auxweights vector can not include any NAs or NaNs or Infs.")
    }
  }
  if (!is.null(signres)){
    if (signres!="pos" & signres!="neg"){
      stop("signres must be either NULL or pos or neg.")
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
  auxdat <- auxdat[,c(pvar,comvar)]

  # add a weight vector to use 'lm' later
  maindat$mainweights <- mainweights
  auxdat$auxweights   <- auxweights

  # number of regressors in a regrssion model (assuming there is only one omitted variable)
  nr <- length(comvar)+1

  #############
  # estimate CDF and Quantile function
  #############

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
  Fypar[is.na(Fypar)] <- 0
  yhat  <- as.matrix(maindat[,comvar])%*%Fypar
  ysd   <- sd(oout1$residuals,na.rm=TRUE)

  # estimate f(pvar | ovar)
  if (ptype==1){

    # continuous proxy variables
    if (is.null(auxweights)){
      pout <- cproxyme(dat=auxdat[,pvar],anchor=1)
    } else{
      pout <- cproxyme(dat=auxdat[,pvar],anchor=1,weights=auxweights)
    }

    if (normalize==TRUE){

      # noramlize proxy variables so that latent variable has mean 0 and std 1
      for (g in 1:length(pvar)){
        auxdat[,pvar[g]] <- (auxdat[,pvar[g]] - pout$mtheta)/(sqrt(pout$vartheta))
      }

      # reestimate measurement equations with normalized proxy variables
      if (is.null(auxweights)){
        pout <- cproxyme(dat=auxdat[,pvar],anchor=1)
      } else{
        pout <- cproxyme(dat=auxdat[,pvar],anchor=1,weights=auxweights)
      }
    }

    alpha0   <- pout$alpha0
    alpha1   <- pout$alpha1
    varnu    <- pout$varnu
    mtheta   <- pout$mtheta
    vartheta <- pout$vartheta

  } else if (ptype==2){
    if (is.null(auxweights)){
      pout <- dproxyme(dat=auxdat[,pvar],sbar,initvar=1)
    } else{
      pout <- dproxyme(dat=auxdat[,pvar],sbar,initvar=1,weights=auxweights)
    }
    M_param     <-pout$M_param
    M_param_col <-pout$M_param_col
    M_param_row <-pout$M_param_row
    mparam      <-pout$mparam
    typeprob    <-pout$typeprob

  } else {
    stop("ptype should be either 1 or 2.")
  }

  N <- dim(auxdat)[1]
  nc <- length(comvar)

  # estimate N(ovar | comvar)
  if (ptype==1){

    # construct normalized proxy variables
    npdat  <- auxdat[,pvar]
    np <- length(pvar)
    nsdnu <- rep(NA,np)

    for (i in 1:np){
      npdat[,i] <- (npdat[,i]-alpha0[i])/alpha1[i]
      nsdnu[i]  <- sqrt(varnu[i]/(alpha1[i]^2))
    }

    # stack up the normalized proxy data
    sdat <- cbind(npdat[,1],auxdat[,comvar])
    colnames(sdat) <- c("y",comvar)
    for (a in 2:np){
      sdat0 <- cbind(npdat[,a],auxdat[,comvar])
      colnames(sdat0) <- c("y",comvar)
      sdat <- rbind(sdat,sdat0)
    }
    sdat <- as.data.frame(sdat)

    f2 <- paste0("y ~ 0 +",comvar[1])
    if (length(comvar)>1){
      for (k in 2:length(comvar)){
        f2 <- paste0(f2,"+",comvar[k])
      }
    }

    if (is.null(auxweights)){
      oout2 <- lm(formula=f2,data=sdat) ## regression without intercept because of "con" in "comvar"
    } else{
      sdat$weights <- rep(auxweights,np)
      oout2 <- lm(formula=f2,data=sdat,weights=weights) ## regression without intercept because of "con" in "comvar"
    }

    # prediction in main data, not auxiliary data
    param <- oout2$coefficients
    param[is.na(param)] <- 0

    Fopar <- matrix(param[1:nc],ncol=1)
    ohat  <- as.matrix(maindat[,comvar])%*%Fopar

    varNoNA <- function(x) var(x,na.rm=TRUE)
    res <- sdat[,"y"] - as.matrix(sdat[,comvar])%*%Fopar
    osd <- mean(sqrt(apply(matrix(res,ncol=np),2,varNoNA)-(nsdnu)^2))

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

  } else if (ptype==2){

    if (is.null(auxweights)){
      oout2 <- multinom(formula=typeprob~as.matrix(auxdat[,comvar[1:(nc-1)]]),maxit=10000,trace=FALSE) ## regression without intercept because of "con" in "comvar"
    } else{
      oout2 <- multinom(formula=typeprob~as.matrix(auxdat[,comvar[1:(nc-1)]]),weights=auxweights,maxit=10000,trace=FALSE) ## regression without intercept because of "con" in "comvar"
    }


    param <- t(coef(oout2))
    param[is.na(param)]<-0

    npr <- dim(param)[1]
    npc <- dim(param)[2]

    # move intercept to the last row
    Fopar <- rbind(matrix(param[2:npr,],ncol=npc),matrix(param[1,],ncol=npc))

    # prediction in main data, not auxiliary data
    Fopar  <- cbind(rep(0,nc),Fopar)
    oprob  <- exp(as.matrix(maindat[,comvar])%*%Fopar)
    oprob  <- oprob/matrix(rep(apply(oprob,1,sum),sbar),ncol=sbar)

    coprob <- t(apply(oprob,1,cumsum))

    #############
    # compute bounds of E[(depvar)*(omitted variable)]
    #############

    ovar_m_l <- rep(NA,Nm)
    ovar_m_u <- rep(NA,Nm)

    if (normalize==TRUE){
      ogrid <- (c(1:sbar)-mean(c(1:sbar)))/sd(c(1:sbar))
    } else{
      ogrid <- c(1:sbar)
    }

    for (k in 1:Nm){
      if (!is.na(maindat[k,depvar]) & !is.nan(maindat[k,depvar]) & !is.na(yhat[k]) & !is.nan(yhat[k]) & !is.na(ysd) & !is.nan(ysd) & sum(is.na(coprob[k,])|is.nan(coprob[k,]))==0 ){
        ovar_m_u[k] <- ogrid[which(   pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd) <coprob[k,])[1]]
        ovar_m_l[k] <- ogrid[which((1-pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd))<coprob[k,])[1]]
      }
    }

  } else {
    stop("ptype should be either 1 or 2.")
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

  if (is.null(mainweights)){

    mu_l <- sum(maindat[,depvar]*ovar_m_l) / sum(Imaindat[,depvar]*Iovar_m_l)
    mu_u <- sum(maindat[,depvar]*ovar_m_u) / sum(Imaindat[,depvar]*Iovar_m_u)

  } else{

    mu_l <- sum(maindat[,depvar]*ovar_m_l*mainweights) / sum(Imaindat[,depvar]*Iovar_m_l*mainweights)
    mu_u <- sum(maindat[,depvar]*ovar_m_u*mainweights) / sum(Imaindat[,depvar]*Iovar_m_u*mainweights)
  }

  # submatrices
  if (ptype==1){

    Inpdat <- !is.na(npdat)
    npdat[!Inpdat] <- 0

    # continuous
    A1 <- vartheta + mtheta^2
    # use normalized proxies to compute covariance, A2
    A2 <- matrix(NA,nrow=1,ncol=nc)
    for (k in 1:nc){
      if (is.null(auxweights)){
        A2[1,k] <- sum(rep(auxdat[,comvar[k]],np)*matrix(as.matrix(npdat),ncol=1)) / sum(rep(Iauxdat[,comvar[k]],np)*matrix(as.matrix(Inpdat),ncol=1))
      } else{
        A2[1,k] <- sum(rep(auxweights*auxdat[,comvar[k]],np)*matrix(as.matrix(npdat),ncol=1)) / sum(rep(auxweights*Iauxdat[,comvar[k]],np)*matrix(as.matrix(Inpdat),ncol=1))
      }
    }
  } else if (ptype==2){

    meanNoNA <- function(x){
      return(mean(x,na.rm=TRUE))
    }

    # discrete
    iprob <- apply(typeprob,2,meanNoNA)
    A1 <- sum((ogrid^2)*iprob)
    A2 <- matrix(0,nrow=1,ncol=nc)
    for (k in 1:nc){
      temp <- 0
      for (l in 1:sbar){
        temp <- temp + ogrid[l]*auxdat[,comvar[k]]*typeprob[,l]
      }
      if (is.null(auxweights)){
        A2[1,k] <- sum(temp) / sum(Iauxdat[,comvar[k]])
      } else{
        A2[1,k] <- sum(temp*auxweights) / sum(Iauxdat[,comvar[k]]*auxweights)
      }
      rm(temp)
    }

  } else{
    stop("ptype must be either 1 or 2")
  }

  if (is.null(auxweights) & is.null(mainweights)){

    C  <- as.matrix(rbind( maindat[,comvar], auxdat[,comvar]))
    IC <- as.matrix(rbind(Imaindat[,comvar],Iauxdat[,comvar]))

    A3 <- (t(C)%*%C)/(t(IC)%*%IC)

  } else if(!is.null(auxweights) & is.null(mainweights)){

    aw <- matrix(rep(auxweights, length(comvar)),ncol=length(comvar)) *(1/sum(auxweights)) * Na

    C   <- as.matrix(rbind( maindat[,comvar],aw* auxdat[,comvar]))
    IC  <- as.matrix(rbind(Imaindat[,comvar],aw*Iauxdat[,comvar]))

    C2  <- as.matrix(rbind( maindat[,comvar], auxdat[,comvar]))
    IC2 <- as.matrix(rbind(Imaindat[,comvar],Iauxdat[,comvar]))

    A3 <- (t(C)%*%C2)/(t(IC)%*%IC2)

  } else if(is.null(auxweights) & !is.null(mainweights)){

    mw <- matrix(rep(mainweights,length(comvar)),ncol=length(comvar)) *(1/sum(mainweights)) * Nm

    C  <- as.matrix(rbind(mw* maindat[,comvar],  auxdat[,comvar]))
    IC <- as.matrix(rbind(mw*Imaindat[,comvar], Iauxdat[,comvar]))

    C2  <- as.matrix(rbind( maindat[,comvar],  auxdat[,comvar]))
    IC2 <- as.matrix(rbind(Imaindat[,comvar], Iauxdat[,comvar]))

    A3 <- (t(C)%*%C2)/(t(IC)%*%IC2)

  } else{

    mw <- matrix(rep(mainweights,length(comvar)),ncol=length(comvar)) *(1/sum(mainweights)) * Nm
    aw <- matrix(rep(auxweights, length(comvar)),ncol=length(comvar)) *(1/sum(auxweights))  * Na

    C  <- as.matrix(rbind(mw* maindat[,comvar], aw* auxdat[,comvar]))
    IC <- as.matrix(rbind(mw*Imaindat[,comvar], aw*Iauxdat[,comvar]))

    C2  <- as.matrix(rbind( maindat[,comvar],  auxdat[,comvar]))
    IC2 <- as.matrix(rbind(Imaindat[,comvar], Iauxdat[,comvar]))

    A3 <- (t(C)%*%C2)/(t(IC)%*%IC2)

  }

  XX <- as.matrix(rbind(cbind(A1,A2),cbind(t(A2),A3)))


  # OLS formula
  if (is.null(mainweights)){
    B <- (t(as.matrix(maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/(t(as.matrix(Imaindat[,depvar]))%*%as.matrix(Imaindat[,comvar]))
  } else {
    B <- (t(as.matrix(mainweights*maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/(t(as.matrix(mainweights*Imaindat[,depvar]))%*%as.matrix(Imaindat[,comvar]))
  }

  B_l <- matrix(c(mu_l,B),ncol=1)
  B_u <- matrix(c(mu_u,B),ncol=1)

  hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
  hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)

  colnames(hat_beta_l) <- c("ovar",comvar)
  colnames(hat_beta_u) <- c("ovar",comvar)

  if (!is.null(signres)){
    if (signres=="pos" & (hat_beta_l[1]<0)){
      # solve the inverse problem
      M <- pinv(XX)
      mu_zero <- -(M[1,2:nr]%*%matrix(B,ncol=1))/M[1,1]

      if (M[1,1]<0){
        mu_u <- mu_zero
        mu_l <- min(mu_zero,mu_l)
      } else{
        mu_l <- mu_zero
        mu_u <- max(mu_zero,mu_u)
      }

      B_l <- matrix(c(mu_l,B),ncol=1)
      B_u <- matrix(c(mu_u,B),ncol=1)

      hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
      hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)

      colnames(hat_beta_l) <- c("ovar",comvar)
      colnames(hat_beta_u) <- c("ovar",comvar)

    }

    if (signres=="neg" & (hat_beta_u[1]>0)){
      # solve the inverse problem
      M <- pinv(XX)
      mu_zero <- -(M[1,2:nr]%*%matrix(B,ncol=1))/M[1,1]

      if (M[1,1]<0){
        mu_l <- mu_zero
        mu_u <- max(mu_zero,mu_u)
      } else{
        mu_u <- mu_zero
        mu_l <- min(mu_zero,mu_l)
      }

      B_l <- matrix(c(mu_l,B),ncol=1)
      B_u <- matrix(c(mu_u,B),ncol=1)

      hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
      hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)

      colnames(hat_beta_l) <- c("ovar",comvar)
      colnames(hat_beta_u) <- c("ovar",comvar)

    }
  }

  # change the order of OLS coefficients
  comvar2 <- comvar[comvar!="con"]
  hat_beta_l <- c(hat_beta_l[,"con"],hat_beta_l[,"ovar"],hat_beta_l[,comvar2])
  hat_beta_u <- c(hat_beta_u[,"con"],hat_beta_u[,"ovar"],hat_beta_u[,comvar2])

  return(list(hat_beta_l=hat_beta_l,hat_beta_u=hat_beta_u,mu_l=mu_l,mu_u=mu_u))
}
