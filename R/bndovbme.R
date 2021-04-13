#' @title bndovbme
#' @description This function runs a two sample least squares when main data contains a dependent variable and 
#' every right hand side regressor but one omitted variable.
#' The function requires an auxiliary data which includes every right hand side regressor but one omitted variable, 
#' and enough proxy variables for the omitted variable. 
#' When the omitted variable is continuous, the auxiliary data must contain at least two continuous proxy variables.
#' When the omitted variable is discrete, the auxiliary data must contain at least three continuous proxy variables.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references Hwang, Yujung (2021). Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.
#' @importFrom utils install.packages
#' @import stats
#' @import np
#' @import pracma
#' @import dplyr
#' @import MASS
#' @import factormodel
#' @import maxLik
#'
#' @param maindat Main data set
#' @param auxdat Auxiliary data set
#' @param depvar A name of a dependent variable in main dataset
#' @param pvar A vector of the names of the proxy variables for the omitted variable. 
#' When proxy variables are continuous, the first proxy variable is used as an anchoring variable.
#' When proxy variables are discrete, the first proxy variable is used for initialization (For details, see a documentation for "dproxyme" function).
#' @param ptype Either 1 (continuous) or 2 (discrete). Whether proxy variables are continuous or discrete. Default is 1 (continuous). 
#' @param comvar A vector of the names of the common regressors existing in both main data and auxiliary data
#' @param sbar A cardinality of the support of the discrete proxy variables. Default is 2. If proxy variables are continuous, this variable is irrelevant.
#' @param coefub An upper bound constraint in the Maximum Likelihood estimation. Default is 100. If you do not want to set any constraint, just set it to a very large number.
#' @param coeflb A lower bound constraint in the Maximum Likelihood estimation coefficients. Default is -100. If you do not want to set any constraint, just set it to a very small number.
#'
#' @return Returns a list of 2 components : \describe{
#' \item{hat_beta_l}{lower bound estimates of regression coefficients}
#'
#' \item{hat_beta_u}{upper bound estimates of regression coefficients}}
#'
#' @export
bndovbme <- function(maindat,auxdat,depvar,pvar,ptype=1,comvar,sbar=2,coefub=100,coeflb=-100){
  
  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("np")
  requireNamespace("pracma")
  requireNamespace("factormodel")
  requireNamespace("maxLik")
  
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
  oout1 <- lm(formula=f1,data=maindat) ## regression without intercept because of "con" in "comvar"
  Fypar <- matrix(oout1$coefficients,ncol=1)
  yhat  <- as.matrix(maindat[,comvar])%*%Fypar
  ysd   <- sd(oout1$residuals,na.rm=TRUE)
    
  # estimate f(pvar | ovar)
  if (ptype==1){
    # continuous proxy variables
    pout <- cproxyme(dat=auxdat[,pvar],anchor=1)
    alpha0   <- pout$alpha0
    alpha1   <- pout$alpha1
    varnu    <- pout$varnu
    mtheta   <- pout$mtheta
    vartheta <- pout$vartheta
    
  } else if (ptype==2){
    pout <- dproxyme(dat=auxdat[,pvar],sbar,initvar=1)
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
    
    Plike <- function(param,cdat,npdat,nsdnu,nc,N){
      
      Fopar <- matrix(param[1:nc],ncol=1)
      osd <- abs(param[(nc+1)]) # positive value
      
      mmu <- as.matrix(cdat)%*%Fopar
      
      ### convolution of two normals
      pdffun <- function(zz,muo,nsdnu) integrate(function(eps,zz,muo,nsdnu) dnorm(zz-eps,mean=muo,sd=osd)*dnorm(eps,mean=0,sd=nsdnu),-Inf,Inf,zz=zz,muo=muo,nsdnu=nsdnu,stop.on.error=FALSE)$value
      
      ll <- rep(0,N)
      for (i in 1:N){
        # convolution
        for (k in 1:np){
          ll[i] <- ll[i] + log(pdffun(zz=npdat[i,k],muo=mmu[i],nsdnu=nsdnu[k]))
        }
        if (is.infinite(ll[i])){
          ll[i] <- -10^(-323) ### lowest number
        }
      }
      return(ll)
    }
    
    # estimate parameters
    A <- rbind( -eye((nc+1)),  eye((nc+1)))
    B <- c(rep(coefub,nc),coefub,-rep(coeflb,nc),-0.001)
    
    oout <- maxLik(logLik=Plike,start=c(rep(0,nc),1),constraints=list(ineqA=A, ineqB=B),cdat=auxdat[,comvar],npdat=npdat,nsdnu=nsdnu,nc=nc,N=N)
    
    
    # prediction in main data, not auxiliary data
    param <- coef(oout)
    Fopar <- matrix(param[1:nc],ncol=1)
    osd   <- abs(param[(nc+1)]) # positive value
    ohat  <- as.matrix(maindat[,comvar])%*%Fopar
    
    #############
    # compute bounds of E[(depvar)*(omitted variable)]
    #############
    
    ovar_m_l <- rep(NA,Nm)
    ovar_m_u <- rep(NA,Nm)
    
    for (k in 1:Nm){
      ovar_m_u[k] <- qnorm(p=   pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd) ,mean=ohat[k],sd=osd)
      ovar_m_l[k] <- qnorm(p=(1-pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd)),mean=ohat[k],sd=osd)
    }
    
  } else if (ptype==2){

    meanNoNA <- function(x){
      return(mean(x,na.rm=TRUE))
    }
    
    # multinomial logit
    Plike <- function(param,cdat,typeprob,sbar,nc,N){
    
      dim(param) <- c(nc,sbar-1)
      param <- cbind(rep(0,nc),param)
      
      nn <- exp(as.matrix(cdat)%*%as.matrix(param))

      dd <- matrix(rep(apply(nn,1,sum),sbar),ncol=sbar)
      pp <- nn/dd
      
      ll <- apply(typeprob*log(pp),1,sum)
      
      return(ll)
    }
    
    # estimate parameters
    A <- rbind( -eye((nc*(sbar-1))),  eye((nc*(sbar-1))))
    B <- c(rep(coefub,nc*(sbar-1)),-rep(coeflb,nc*(sbar-1)))
    
    oout <- maxLik(logLik=Plike,start=c(rep(0,nc*(sbar-1))),constraints=list(ineqA=A, ineqB=B),cdat=auxdat[,comvar],typeprob=typeprob,sbar=sbar,nc=nc,N=N)
    
    # prediction in main data, not auxiliary data
    Fopar  <- matrix(coef(oout),ncol=(sbar-1))
    Fopar  <- cbind(rep(0,nc),Fopar)
    oprob  <- exp(as.matrix(maindat[,comvar])%*%Fopar)
    oprob  <- oprob/matrix(rep(apply(oprob,1,sum),sbar),ncol=sbar)
    
    coprob <- t(apply(oprob,1,cumsum))
    
    #############
    # compute bounds of E[(depvar)*(omitted variable)]
    #############
    
    ovar_m_l <- rep(NA,Nm)
    ovar_m_u <- rep(NA,Nm)
    
    for (k in 1:Nm){
      ovar_m_u[k] <- which(   pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd) <coprob[k,])[1]
      ovar_m_l[k] <- which((1-pnorm(q=maindat[k,depvar],mean=yhat[k],sd=ysd))<coprob[k,])[1]
    }
    
  } else {
    stop("ptype should be either 1 or 2.")    
  }
  
  #############
  # compute lower bound and upper bound
  #############
  
  mu_l <- mean(maindat[,depvar]*ovar_m_l,na.rm=TRUE)
  mu_u <- mean(maindat[,depvar]*ovar_m_u,na.rm=TRUE)
  
  hat_beta_l <- rep(NA,nr)
  hat_beta_u <- rep(NA,nr)
  
  # submatrices
  if (ptype==1){
    # continuous
    A1 <- vartheta + mtheta^2
    # use normalized proxies to compute covariance, A2
    A2 <- matrix(NA,nrow=1,ncol=nc)
    for (k in 1:nc){
      A2[1,k] <- mean(rep(auxdat[,comvar[k]],np)*matrix(as.matrix(npdat),ncol=1), na.rm=TRUE) 
    }
  } else if (ptype==2){
    
    # discrete
    iprob <- apply(typeprob,2,meanNoNA)
    A1 <- sum((c(1:sbar)^2)*iprob)
    A2 <- matrix(0,nrow=1,ncol=nc)
    for (k in 1:nc){
      temp <- 0
      for (l in 1:sbar){
        temp <- temp + l*auxdat[,comvar[k]]*typeprob[,l]        
      }
      A2[1,k] <- mean(temp,na.rm=TRUE)
      rm(temp)
    }
    
  } else{
    stop("ptype must be either 1 or 2")
  }
  
  C  <- as.matrix(rbind(maindat[,comvar],auxdat[,comvar]))
  A3 <- (t(C)%*%C)/(Na+Nm)
  
  XX <- as.matrix(rbind(cbind(A1,A2),cbind(t(A2),A3)))
  
  # OLS formula
  B <- (t(as.matrix(maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/Nm
  B_l <- matrix(c(mu_l,B),ncol=1)
  B_u <- matrix(c(mu_u,B),ncol=1)
  
  hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
  hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
  
  colnames(hat_beta_l) <- c("ovar",comvar)
  colnames(hat_beta_u) <- c("ovar",comvar)
  
  # change the order of OLS coefficients
  comvar2 <- comvar[comvar!="con"]
  hat_beta_l <- c(hat_beta_l[,"con"],hat_beta_l[,"ovar"],hat_beta_l[,comvar2])
  hat_beta_u <- c(hat_beta_u[,"con"],hat_beta_u[,"ovar"],hat_beta_u[,comvar2])
  
  return(list(hat_beta_l=hat_beta_l,hat_beta_u=hat_beta_u))
}
