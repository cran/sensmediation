#'Example data for the functions in sensmediation
#'
#'The data are a subsample of 1000 observations from Riksstroke, the Swedish Stroke Register. The original data consisted
#'of over 50 000 patients with first time ischemic stroke during the years 2009-2012. The data are limited to patients over
#'the age of 44 and its purpose is to illustrate the functioning of the functions in the package.
#'
#'@docType data
#'@keywords datasets
#'
#'@usage data(RSdata)
#'
#'@format A data frame with 1000 observations on the following 5 variables.
#' \describe{
#'   \item{\code{cf.3mo}}{Outcome: case fatality within 3 months after stroke, 1 = deceased, 0 = not deceased.}
#'   \item{\code{lowered.consc}}{Mediator: level of consciousness upon arrival to hospital. 1 = lowered consciousness, 0 = fully alert.}
#'   \item{\code{AF}}{Exposure: atrial fibrillation. Factor with levels, "1" = atrial fibrillation, "0" = no atrial fibrillation.}
#'   \item{\code{age.cat}}{Age at time of stroke. Factor with levels, "45-69", "70-79", "80-89" and "90-".}
#'   \item{\code{sex}}{Factor with levels, "1" = male, "0" = female}
#' }
#' 
#' 
#' @examples 
#' data(RSdata)
"RSdata"


