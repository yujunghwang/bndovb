#' @title bndovbme_tuning
#' @description This function computes an optimal tuning parameter to compute the confidence interval for bndovbme function
#' The function returns an optimal tuning parameter using double bootstrap procedure
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Available at SSRN.\doi{10.2139/ssrn.3866876}}}
#' @importFrom utils install.packages
#' @import stats
#' @import np
#' @importFrom pracma pinv eye randi
#' @importFrom MASS mvrnorm
#' @import factormodel
#' @importFrom nnet multinom
#' @import doParallel
#' @import foreach
#' @import parallel
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
#' @param nboot Number of bootstraps to compute the confidence interval. Default is 100.
#' @param scalegrid Tuning parameter grid to search. It must be a vector of numbers between -1/2 and 0. Default is c(-1/2,-1/3,-1/4,-1/5,-1/6).
#' @param tau Significance level. (1-tau)% confidence interval is computed. Default is 0.05.
#' @param seed Seed for random number generation. Default is 210823.
#' @param parallel Either TRUE or FALSE. Whether to compute in parallel. Default is TRUE.
#'
#' @return Returns a list of 3 components : \describe{
#'
#' \item{optimal_scale}{An optimal scale parameter which gives coverage rates closest to (1-tau)}
#'
#' \item{cover_beta_l}{A matrix of coverage rates of the lower bound parameters under different scale parameters}
#'
#' \item{cover_beta_u}{A matrix of coverage rates of the lower bound parameters under different scale parameters}}
#'
#' @examples
#' ## load example data
#' data(maindat_mecont)
#' data(auxdat_mecont)
#'
#' ## set ptype=1 for continuous proxy variables
#' pvar<-c("z1","z2","z3")
#' cvar<-c("x","w1")
#'
#' # To shorten computation time, I set the number of bootstrap small in an example below.
#' # In practice, please set it a large number
#' bndovbme_tuning(maindat_mecont,auxdat_mecont,depvar="y",pvar=pvar,ptype=1,comvar=cvar,nboot=2)
#'
#' @export
bndovbme_tuning <- function(maindat,auxdat,depvar,pvar,ptype=1,comvar,sbar=2,mainweights=NULL,auxweights=NULL,normalize=TRUE,signres=NULL,nboot=100,scalegrid=c(-1/2,-1/3,-1/4,-1/5,-1/6),tau=0.05,seed=210823,parallel=TRUE){

  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("np")
  requireNamespace("pracma")
  requireNamespace("factormodel")
  requireNamespace("nnet")
  requireNamespace("doParallel")
  requireNamespace("foreach")
  requireNamespace("parallel")

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

  if (nboot<2){
    stop("The number of bootstrap is too small. Enter a number greater than 1.")
  }

  if ( sum((scalegrid < -1/2) | (scalegrid > 0)) ){
    stop("The scale parameters in the scalegrid must be between -1/2 and 0.")
  }

  if ((tau<0) | (tau>1)){
    stop("tau must be between 0 and 1.")
  }

  ###############

  # estimate once
  oout <- bndovbme(maindat=maindat,auxdat=auxdat,depvar=depvar,pvar=pvar,ptype=ptype,comvar=comvar,sbar=sbar,mainweights=mainweights,auxweights=auxweights,normalize=normalize,signres=signres,ci=FALSE,display=FALSE)

  hat_beta_l <- oout[[1]]
  hat_beta_u <- oout[[2]]

  ne <- length(scalegrid)

  cover_beta_l <- array(NA,dim=c(length(hat_beta_l),ne,nboot))
  cover_beta_u <- array(NA,dim=c(length(hat_beta_u),ne,nboot))

  #############
  # Double bootstrap
  #############

  set.seed(seed)

  # number of observations
  Nm <- dim(maindat)[1]
  Na <- dim(auxdat)[1]

  # draw bootstrap samples with replacement
  bmain_ind <- randi(Nm,n=Nm,m=nboot)
  baux_ind  <- randi(Na,n=Na,m=nboot)

  # progress message
  prog <- round(quantile(c(1:nboot),probs=seq(0.1,1,0.1)),digits=0)

  if (parallel==FALSE){
    # compute serially
    for (b1 in 1:nboot){

      if (b1%in%prog){
        print(paste0(names(prog)[max(which(prog==b1))]," completed"))
      }

      # bootstrap sample
      bmaindat <- maindat[bmain_ind[,b1],]
      bauxdat  <- auxdat[ baux_ind[,b1],]

      if (!is.null(mainweights)){
        bmainweights <- mainweights[bmain_ind[,b1]]
      } else{
        bmainweights <- NULL
      }

      if (!is.null(auxweights)){
        bauxweights <- auxweights[baux_ind[,b1]]
      } else{
        bauxweights <- NULL
      }

      for (e1 in 1:ne){

        # compute confidence intervals by drawing second bootstrap. Use different seed number.
        boout <- bndovbme(maindat=bmaindat,auxdat=bauxdat,depvar=depvar,pvar=pvar,ptype=ptype,comvar=comvar,sbar=sbar,
                          mainweights=bmainweights,auxweights=bauxweights,normalize=normalize,signres=signres,
                          ci=TRUE,nboot=nboot,scale=scalegrid[e1],tau=tau,seed=b1,display=FALSE)

        # confidence intervals
        bhat_beta_l_cil <- boout[[5]]
        bhat_beta_l_ciu <- boout[[6]]
        bhat_beta_u_cil <- boout[[7]]
        bhat_beta_u_ciu <- boout[[8]]

        # count the coverage
        cover_beta_l[,e1,b1] <- as.numeric((bhat_beta_l_cil<=hat_beta_l) & (bhat_beta_l_ciu>=hat_beta_l))
        cover_beta_u[,e1,b1] <- as.numeric((bhat_beta_u_cil<=hat_beta_u) & (bhat_beta_u_ciu>=hat_beta_u))
      }
    }
  } else{

    # compute in parallel
    print("bootstrap in progress in parallel...")
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      ncore <- 2L
    } else {
      # use all cores in devtools::test()
      ncore <- parallel::detectCores()-1
    }
    print(paste0("Using ",ncore," workers..."))
    cl<-makeCluster(ncore, setup_strategy = "sequential")
    registerDoParallel(cl)
    fout <- foreach(b1=1:nboot,.packages = c("stats","utils","np","pracma","factormodel","nnet"),.export=ls(globalenv())) %dopar% {

      # bootstrap sample
      bmaindat <- maindat[bmain_ind[,b1],]
      bauxdat  <- auxdat[ baux_ind[,b1],]

      if (!is.null(mainweights)){
        bmainweights <- mainweights[bmain_ind[,b1]]
      } else{
        bmainweights <- NULL
      }

      if (!is.null(auxweights)){
        bauxweights <- auxweights[baux_ind[,b1]]
      } else{
        bauxweights <- NULL
      }

      temp1 <- array(NA,dim=c(length(hat_beta_l),ne))
      temp2 <- array(NA,dim=c(length(hat_beta_u),ne))

      for (e1 in 1:ne){
        # compute confidence intervals by drawing second bootstrap. Use different seed number.
        boout <- bndovbme(maindat=bmaindat,auxdat=bauxdat,depvar=depvar,pvar=pvar,ptype=ptype,comvar=comvar,sbar=sbar,
                          mainweights=bmainweights,auxweights=bauxweights,normalize=normalize,signres=signres,
                          ci=TRUE,nboot=nboot,scale=scalegrid[e1],tau=tau,seed=b1,display=FALSE)

        # confidence intervals
        bhat_beta_l_cil <- boout[[5]]
        bhat_beta_l_ciu <- boout[[6]]
        bhat_beta_u_cil <- boout[[7]]
        bhat_beta_u_ciu <- boout[[8]]

        # count the coverage
        temp1[,e1] <- as.numeric((bhat_beta_l_cil<=hat_beta_l) & (bhat_beta_l_ciu>=hat_beta_l))
        temp2[,e1] <- as.numeric((bhat_beta_u_cil<=hat_beta_u) & (bhat_beta_u_ciu>=hat_beta_u))
      }
      temp <- rbind(temp1,temp2)

      return(temp)
    }
    stopCluster(cl)
    fout2 <- array(unlist(fout),dim=c((length(hat_beta_l)+length(hat_beta_u)),ne,nboot))
    cover_beta_l <- fout2[1:length(hat_beta_l),,]
    cover_beta_u <- fout2[(length(hat_beta_l)+1):(length(hat_beta_l)+length(hat_beta_u)),,]
  }

  # average the coverage rate
  cover_beta_l <- apply(cover_beta_l,c(1,2),mean)
  cover_beta_u <- apply(cover_beta_u,c(1,2),mean)

  colnames(cover_beta_l) <- scalegrid
  colnames(cover_beta_u) <- scalegrid
  rownames(cover_beta_l) <- names(hat_beta_l)
  rownames(cover_beta_u) <- names(hat_beta_u)

  aa <- rbind(cover_beta_l,cover_beta_u)
  acover <- apply(aa,2,mean)

  dif <- abs(acover - (1-tau))

  # choose the optimal scale parameter
  optimal_scale_ind <- which(dif==min(dif))[1]
  optimal_scale     <- scalegrid[optimal_scale_ind]

  return(list(optimal_scale=optimal_scale,cover_beta_l=cover_beta_l,cover_beta_u=cover_beta_u))
}
