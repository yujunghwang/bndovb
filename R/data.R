#' A simulated main data to show how to use 'bndovb' function
#'
#' @format A data frame with 100000 rows and 3 variables:
#' \describe{
#'   \item{x2}{A common covariate in both main and auxiliary data}
#'   \item{x3}{A common covariate in both main and auxiliary data}
#'   \item{y}{A dependent variable}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"maindat_nome"

#' A simulated auxiliary data to show how to use 'bndovb' function
#'
#' @format A data frame with 50000 rows and 3 variables:
#' \describe{
#'   \item{x1}{An omitted variable in the main data}
#'   \item{x2}{A common covariate in both main and auxiliary data}
#'   \item{x3}{A common covariate in both main and auxiliary data}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"auxdat_nome"

#' A simulated main data to show how to use 'bndovbme' function with continuous proxy variables
#'
#' @format A data frame with 3000 rows and 3 variables:
#' \describe{
#'   \item{w1}{A common covariate in both main and auxiliary data}
#'   \item{x}{A common covariate in both main and auxiliary data}
#'   \item{y}{A dependent variable}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"maindat_mecont"

#' A simulated auxiliary data to show how to use 'bndovbme' function with continuous proxy variables
#'
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'   \item{w1}{A common covariate in both main and auxiliary data}
#'   \item{x}{A common covariate in both main and auxiliary data}
#'   \item{z1}{A continuous proxy variable}
#'   \item{z2}{A continuous proxy variable}
#'   \item{z3}{A continuous proxy variable}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"auxdat_mecont"

#' A simulated main data to show how to use 'bndovbme' function with discrete proxy variables
#'
#' @format A data frame with 3000 rows and 3 variables:
#' \describe{
#'   \item{w1}{A common covariate in both main and auxiliary data}
#'   \item{x}{A common covariate in both main and auxiliary data}
#'   \item{y}{A dependent variable}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"maindat_medisc"

#' A simulated auxiliary data to show how to use 'bndovbme' function with discrete proxy variables
#'
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'   \item{w1}{A common covariate in both main and auxiliary data}
#'   \item{x}{A common covariate in both main and auxiliary data}
#'   \item{z1}{A discrete proxy variable}
#'   \item{z2}{A discrete proxy variable}
#'   \item{z3}{A discrete proxy variable}
#' }
#' @source This dataset was simulated by simulatePackageData.R in data-raw folder
"auxdat_medisc"


