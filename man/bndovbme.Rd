% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bndovbme.R
\name{bndovbme}
\alias{bndovbme}
\title{bndovbme}
\usage{
bndovbme(
  maindat,
  auxdat,
  depvar,
  pvar,
  ptype = 1,
  comvar,
  sbar = 2,
  mainweights = NULL,
  auxweights = NULL,
  normalize = TRUE,
  signres = NULL,
  ci = FALSE,
  nboot = 100,
  scale = -1/2,
  tau = 0.05,
  seed = 210823,
  display = TRUE
)
}
\arguments{
\item{maindat}{Main data set. It must be a data frame.}

\item{auxdat}{Auxiliary data set. It must be a data frame.}

\item{depvar}{A name of a dependent variable in main dataset}

\item{pvar}{A vector of the names of the proxy variables for the omitted variable.
When proxy variables are continuous, the first proxy variable is used as an anchoring variable.
When proxy variables are discrete, the first proxy variable is used for initialization (For details, see a documentation for "dproxyme" function).}

\item{ptype}{Either 1 (continuous) or 2 (discrete). Whether proxy variables are continuous or discrete. Default is 1 (continuous).}

\item{comvar}{A vector of the names of the common regressors existing in both main data and auxiliary data}

\item{sbar}{A cardinality of the support of the discrete proxy variables. Default is 2. If proxy variables are continuous, this variable is irrelevant.}

\item{mainweights}{An optional weight vector for the main dataset. The length must be equal to the number of rows of 'maindat'.}

\item{auxweights}{An optional weight vector for the auxiliary dataset. The length must be equal to the number of rows of 'auxdat'.}

\item{normalize}{Whether to normalize the omitted variable to have mean 0 and standard deviation 1. Set TRUE or FALSE.
Default is TRUE. If FALSE, then the scale of the omitted variable is anchored with the first proxy variable in pvar list.}

\item{signres}{An option to impose a sign restriction on a coefficient of an omitted variable. Set either NULL or pos or neg.
Default is NULL. If NULL, there is no sign restriction.
If 'pos', the estimator imposes an extra restriction that the coefficient of an omitted variable must be positive.
If 'neg', the estimator imposes an extra restriction that the coefficient of an omitted variable must be negative.}

\item{ci}{An option to compute an equal-tailed confidence interval. Default is FALSE. It may take some time to compute CI from bootstrap.}

\item{nboot}{Number of bootstraps to compute the confidence interval. Default is 100.}

\item{scale}{A tuning parameter for rescaled numerical bootstrap. The value must be between -1/2 and 0. (main data sample size)^scale is the tuning parameter epsilon_n in Hwang (2021). Default is -1/2 (that is, standard bootstrap).}

\item{tau}{Significance level. (1-tau)\% confidence interval is computed. Default is 0.05.}

\item{seed}{Seed for random number generation. Default is 210823.}

\item{display}{It must be either TRUE or FALSE. Whether to display progress and messages. Default is TRUE.}
}
\value{
Returns a list of 4 components : \describe{
\item{hat_beta_l}{lower bound estimates of regression coefficients}

\item{hat_beta_u}{upper bound estimates of regression coefficients}

\item{mu_l}{lower bound estimate of E[ovar*depvar]}

\item{mu_u}{upper bound estimate of E[ovar*depvar]}

\item{hat_beta_l_cil}{(1-tau)\% confidence interval lower bound for hat_beta_l}

\item{hat_beta_l_ciu}{(1-tau)\% confidence interval upper bound for hat_beta_l}

\item{hat_beta_u_cil}{(1-tau)\% confidence interval lower bound for hat_beta_u}

\item{hat_beta_u_ciu}{(1-tau)\% confidence interval upper bound for hat_beta_u}

\item{mu_l_cil}{(1-tau)\% confidence interval lower bound for mu_l}

\item{mu_l_ciu}{(1-tau)\% confidence interval upper bound for mu_l}

\item{mu_u_cil}{(1-tau)\% confidence interval lower bound for mu_u}

\item{mu_u_ciu}{(1-tau)\% confidence interval upper bound for mu_u}}
}
\description{
This function runs a two sample least squares when main data contains a dependent variable and
every right hand side regressor but one omitted variable.
The function requires an auxiliary data which includes every right hand side regressor but one omitted variable,
and enough proxy variables for the omitted variable.
When the omitted variable is continuous, the auxiliary data must contain at least two continuous proxy variables.
When the omitted variable is discrete, the auxiliary data must contain at least three continuous proxy variables.
}
\examples{
## load example data
data(maindat_mecont)
data(auxdat_mecont)

## set ptype=1 for continuous proxy variables
 pvar<-c("z1","z2","z3")
 cvar<-c("x","w1")
bndovbme(maindat=maindat_mecont,auxdat=auxdat_mecont,depvar="y",pvar=pvar,ptype=1,comvar=cvar)

## set ptype=2 for discrete proxy variables
data(maindat_medisc)
data(auxdat_medisc)
bndovbme(maindat=maindat_medisc,auxdat=auxdat_medisc,depvar="y",pvar=pvar,ptype=2,comvar=cvar)

}
\references{
\describe{
\item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Available at SSRN. \doi{10.2139/ssrn.3866876}}}
}
\author{
Yujung Hwang, \email{yujungghwang@gmail.com}
}
