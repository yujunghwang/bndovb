#' @title bndovb
#' @description This function runs a two sample least squares when auxiliary data contains every right-hand side regressor
#' and main data contains a dependent variable and every right-hand side regressor but one omitted variable.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Available at SSRN.\doi{10.2139/ssrn.3866876}}}
#' @importFrom utils install.packages
#' @import stats
#' @import np
#' @importFrom pracma pinv eye randi
#' @importFrom MASS mvrnorm
#'
#' @param maindat Main data set. It must be a data frame.
#' @param auxdat Auxiliary data set. It must be a data frame.
#' @param depvar A name of a dependent variable in main dataset
#' @param ovar A name of an omitted variable in main dataset which exists in auxiliary data
#' @param comvar A vector of the names of common regressors existing in both main data and auxiliary data
#' @param method CDF and Quantile function estimation method.
#' Users can choose either 1 or 2. If the method is 1, the CDF and quantile function is estimated assuming a parametric normal distribution.
#' If the method is 2, the CDF and quantile function is estimated using a nonparaemtric estimator in Li and Racine(2008) \doi{10.1198/073500107000000250}, Li, Lin, and Racine(2013) \doi{10.1080/07350015.2012.738955}.
#' Default is 1.
#' @param mainweights An optional weight vector for the main dataset. The length must be equal to the number of rows of 'maindat'.
#' @param auxweights An optional weight vector for the auxiliary dataset. The length must be equal to the number of rows of 'auxdat'.
#' @param signres An option to impose a sign restriction on a coefficient of an omitted variable. Set either NULL or pos or neg.
#' Default is NULL. If NULL, there is no sign restriction.
#' If 'pos', the estimator imposes an extra restriction that the coefficient of an omitted variable must be positive.
#' If 'neg', the estimator imposes an extra restriction that the coefficient of an omitted variable must be negative.
#' @param ci An option to compute an equal-tailed confidence interval. Default is FALSE. It may take some time to compute CI from bootstrap.
#' @param nboot Number of bootstraps to compute the confidence interval. Default is 100.
#' @param scale A tuning parameter for rescaled numerical bootstrap. The value must be between -1/2 and 0. (main data sample size)^scale is the tuning parameter epsilon_n in Hwang (2021). Default is -1/2 (that is, standard bootstrap).
#' @param tau Significance level. (1-tau)% confidence interval is computed. Default is 0.05.
#' @param seed Seed for random number generation. Default is 210823.
#' @param display It must be either TRUE or FALSE. Whether to display progress and messages. Default is TRUE.
#'
#' @return Returns a list of 12 components : \describe{
#' \item{hat_beta_l}{lower bound estimates of regression coefficients}
#'
#' \item{hat_beta_u}{upper bound estimates of regression coefficients}
#'
#' \item{mu_l}{lower bound estimate of E\[ovar*depvar\]}
#'
#' \item{mu_u}{upper bound estimate of E\[ovar*depvar\]}
#'
#' \item{hat_beta_l_cil}{(1-tau)% confidence interval lower bound for hat_beta_l}
#'
#' \item{hat_beta_l_ciu}{(1-tau)% confidence interval upper bound for hat_beta_l}
#'
#' \item{hat_beta_u_cil}{(1-tau)% confidence interval lower bound for hat_beta_u}
#'
#' \item{hat_beta_u_ciu}{(1-tau)% confidence interval upper bound for hat_beta_u}
#'
#' \item{mu_l_cil}{(1-tau)% confidence interval lower bound for mu_l}
#'
#' \item{mu_l_ciu}{(1-tau)% confidence interval upper bound for mu_l}
#'
#' \item{mu_u_cil}{(1-tau)% confidence interval lower bound for mu_u}
#'
#' \item{mu_u_ciu}{(1-tau)% confidence interval upper bound for mu_u}}
#'
#' @examples
#' data(maindat_nome)
#' data(auxdat_nome)
#'
#' bndovb(maindat=maindat_nome,auxdat=auxdat_nome,depvar="y",ovar="x1",comvar=c("x2","x3"),method=1)
#'
#'
#' @export
bndovb <- function(maindat,auxdat,depvar,ovar,comvar,method=1,mainweights=NULL,auxweights=NULL,signres=NULL,ci=FALSE,nboot=100,scale=-1/2,tau=0.05,seed=210823,display=TRUE){

  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("np")
  requireNamespace("pracma")

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

  if (nboot<2){
    stop("The number of bootstrap is too small. Enter a number greater than 1.")
  }

  if ((scale < -1/2) | (scale > 0)){
    stop("The scale parameter must be between -1/2 and 0.")
  }

  if ((tau<0) | (tau>1)){
    stop("tau must be between 0 and 1.")
  }

  if (!is.logical(ci)){
    stop("ci must be either TRUE or FALSE.")
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

  # add a weight vector to use 'lm' later
  maindat$mainweights <- mainweights
  auxdat$auxweights   <- auxweights

  # number of regressors in a regrssion model
  nr <- length(comvar)+length(ovar)

  # a subroutine computing XX, B_l, B_u, mu_l, mu_u
  bndovb_moments <- function(maindat,auxdat,mainweights,auxweights){

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
      Fypar[is.na(Fypar)] <- 0

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
      Fopar[is.na(Fopar)] <-0

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

      A1 <- (t(as.matrix(auxweights*auxdat[,ovar]))%*%as.matrix(auxdat[,ovar]))  /t(as.matrix(auxweights*Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,ovar])
      A2 <- (t(as.matrix(auxweights*auxdat[,ovar]))%*%as.matrix(auxdat[,comvar]))/t(as.matrix(auxweights*Iauxdat[,ovar]))%*%as.matrix(Iauxdat[,comvar])

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
    } else{
      B <- (t(as.matrix(mainweights*maindat[,depvar]))%*%as.matrix(maindat[,comvar]))/(t(as.matrix(mainweights*Imaindat[,depvar]))%*%as.matrix(Imaindat[,comvar]))
    }

    B_l <- matrix(c(mu_l,B),ncol=1)
    B_u <- matrix(c(mu_u,B),ncol=1)

    return(list(XX=XX,B_l=B_l,B_u=B_u,mu_l=mu_l,mu_u=mu_u))

  }

  # compute XX, B_l, B_u
  mout <- bndovb_moments(maindat,auxdat,mainweights,auxweights)
  XX   <- mout[[1]]
  B_l  <- mout[[2]]
  B_u  <- mout[[3]]
  mu_l <- mout[[4]]
  mu_u <- mout[[5]]

  # subroutine to compute hat_beta_l and hat_beta_u and mu_l and mu_u (sign restriction adjustment) given XX, B_l, B_u, mu_l, mu_u
  # return hat_beta_l, hat_beta_u
  bndovb_coef <- function(XX,B_l,B_u,mu_l,mu_u){

    hat_beta_l <- matrix(pmin(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)
    hat_beta_u <- matrix(pmax(pinv(XX)%*%B_l,pinv(XX)%*%B_u),nrow=1)

    colnames(hat_beta_l) <- c(ovar,comvar)
    colnames(hat_beta_u) <- c(ovar,comvar)

    if (!is.null(signres)){

      # length(ovar)=1
      B <- B_l[2:nr]

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

        # sign restricted model
        rB_l <- matrix(c(mu_l,B),ncol=1)
        rB_u <- matrix(c(mu_u,B),ncol=1)

        hat_beta_l <- matrix(pmin(pinv(XX)%*%rB_l,pinv(XX)%*%rB_u),nrow=1)
        hat_beta_u <- matrix(pmax(pinv(XX)%*%rB_l,pinv(XX)%*%rB_u),nrow=1)

        colnames(hat_beta_l) <- c(ovar,comvar)
        colnames(hat_beta_u) <- c(ovar,comvar)

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

        # sign restricted model
        rB_l <- matrix(c(mu_l,B),ncol=1)
        rB_u <- matrix(c(mu_u,B),ncol=1)

        hat_beta_l <- matrix(pmin(pinv(XX)%*%rB_l,pinv(XX)%*%rB_u),nrow=1)
        hat_beta_u <- matrix(pmax(pinv(XX)%*%rB_l,pinv(XX)%*%rB_u),nrow=1)

        colnames(hat_beta_l) <- c(ovar,comvar)
        colnames(hat_beta_u) <- c(ovar,comvar)

      }
    }

    # change the order of OLS coefficients
    comvar2 <- comvar[comvar!="con"]
    hat_beta_l <- c(hat_beta_l[,"con"],hat_beta_l[,ovar],hat_beta_l[,comvar2])
    hat_beta_u <- c(hat_beta_u[,"con"],hat_beta_u[,ovar],hat_beta_u[,comvar2])

    return(list(hat_beta_l,hat_beta_u,mu_l,mu_u))
  }

  moout2 <- bndovb_coef(XX,B_l,B_u,mu_l,mu_u)

  hat_beta_l <- moout2[[1]]
  hat_beta_u <- moout2[[2]]
  mu_l       <- moout2[[3]]
  mu_u       <- moout2[[4]]

  ###################################
  # Confidence Interval computation
  ###################################

  hat_beta_l_cil <- NULL
  hat_beta_l_ciu <- NULL
  hat_beta_u_cil <- NULL
  hat_beta_u_ciu <- NULL

  mu_l_cil <- NULL
  mu_l_ciu <- NULL
  mu_u_cil <- NULL
  mu_u_ciu <- NULL

  if (ci==TRUE){

    # set seed
    set.seed(seed)

    # draw bootstrap samples with replacement
    bmain_ind <- randi(Nm,n=Nm,m=nboot)
    baux_ind  <- randi(Na,n=Na,m=nboot)

    # matrices to save derivatives
    dhat_beta_l <- array(NA,dim=c(nr,nboot))
    dhat_beta_u <- array(NA,dim=c(nr,nboot))
    dmu_l       <- rep(NA,nboot)
    dmu_u       <- rep(NA,nboot)

    # progress message
    prog <- round(quantile(c(1:nboot),probs=seq(0.1,1,0.1)),digits=0)
    prog_ind <- 1

    for (b1 in 1:nboot){

      if (display==TRUE){
        if (b1%in%prog){
          print(paste0(names(prog)[prog_ind]," completed"))
          prog_ind <- prog_ind + 1
        }
      }

      # bootstrap sample
      bmaindat <- maindat[bmain_ind[,b1],]
      bauxdat  <- auxdat[ baux_ind[,b1],]

      # compute bootstrap moments (return  : XX, B_l, B_u)
      bmout <- bndovb_moments(maindat=bmaindat,auxdat=bauxdat,mainweights=as.vector(bmaindat$mainweights),auxweights=as.vector(bauxdat$auxweights))
      bXX   <- bmout[[1]]
      bB_l  <- bmout[[2]]
      bB_u  <- bmout[[3]]
      bmu_l <- bmout[[4]]
      bmu_u <- bmout[[5]]

      # rescale by sample size

      # tuning parameter
      en <- Nm^scale
      rn <- sqrt(Nm)

      adjbXX   <- XX   + en*rn*(bXX-XX)
      adjbB_l  <- B_l  + en*rn*(bB_l-B_l)
      adjbB_u  <- B_u  + en*rn*(bB_u-B_u)
      adjbmu_l <- mu_l + en*rn*(bmu_l-mu_l)
      adjbmu_u <- mu_u + en*rn*(bmu_u-mu_u)

      # take the derivative
      bmoout2 <- bndovb_coef(adjbXX,adjbB_l,adjbB_u,adjbmu_l,adjbmu_u)

      bhat_beta_l <- bmoout2[[1]]
      bhat_beta_u <- bmoout2[[2]]
      bmu_l       <- bmoout2[[3]]
      bmu_u       <- bmoout2[[4]]

      dhat_beta_l[,b1] <- (bhat_beta_l - hat_beta_l)/en
      dhat_beta_u[,b1] <- (bhat_beta_u - hat_beta_u)/en
      dmu_l[b1]        <- (bmu_l       - mu_l)/en
      dmu_u[b1]        <- (bmu_u       - mu_u)/en

    }

    rquantile <- function(x){
      return(quantile(x,probs=(1-tau/2)))
    }
    lquantile <- function(x){
      return(quantile(x,probs=(tau/2)))
    }

    # find the tau and (1-tau) percentile
    dhat_beta_l_r <- apply(dhat_beta_l,1,rquantile)
    dhat_beta_l_l <- apply(dhat_beta_l,1,lquantile)
    dhat_beta_u_r <- apply(dhat_beta_u,1,rquantile)
    dhat_beta_u_l <- apply(dhat_beta_u,1,lquantile)
    dmu_l_r       <- rquantile(dmu_l)
    dmu_l_l       <- lquantile(dmu_l)
    dmu_u_r       <- rquantile(dmu_u)
    dmu_u_l       <- lquantile(dmu_u)

    # compute the bound
    hat_beta_l_cil <- hat_beta_l - dhat_beta_l_r / rn
    hat_beta_l_ciu <- hat_beta_l - dhat_beta_l_l / rn
    hat_beta_u_cil <- hat_beta_u - dhat_beta_u_r / rn
    hat_beta_u_ciu <- hat_beta_u - dhat_beta_u_l / rn

    mu_l_cil <- mu_l - dmu_l_r /rn
    mu_l_ciu <- mu_l - dmu_l_l /rn
    mu_u_cil <- mu_u - dmu_u_r /rn
    mu_u_ciu <- mu_u - dmu_u_l /rn

  }

  if ((ci==FALSE) & (display==TRUE)){
    print("If you want to compute an equal-tailed confidence interval using a numerical delta method, set ci=TRUE instead. Default is 95% CI. If you want a different coverage, set a different tau.")
  }

  return(list(hat_beta_l=hat_beta_l,hat_beta_u=hat_beta_u,mu_l=mu_l,mu_u=mu_u,
              hat_beta_l_cil=hat_beta_l_cil,hat_beta_l_ciu=hat_beta_l_ciu,hat_beta_u_cil=hat_beta_u_cil,hat_beta_u_ciu=hat_beta_u_ciu,
              mu_l_cil=mu_l_cil,mu_l_ciu=mu_l_ciu,mu_u_cil=mu_u_cil,mu_u_ciu=mu_u_ciu))
}
