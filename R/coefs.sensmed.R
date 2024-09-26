#'ML estimation of regression parameters for calculation of direct and indirect effects under unobserved confounding
#'
#'This function gives ML estimates of the regression parameters used to calculate mediation effects and perform sensitivity analysis. The optimization is
#'performed using \code{\link[maxLik]{maxLik}}, see Details for more information. Called by \code{\link{sensmediation}}.
#'@param model.expl Fitted \code{\link{glm}} model object. If sensitivity analysis to mediator-outcome confounding the mediator model. Otherwise the exposure model.
#'@param model.resp Fitted \code{\link{glm}} model object. If sensitivity analysis to exposure-mediator confounding the mediator model. Otherwise the outcome model.
#'@param Rho The sensitivity parameter vector. If \code{type="my"} the correlation between the error terms in the mediator and outcome models. If \code{type="zm"} the correlation between the error terms in the exposure and mediator models. If \code{type="zy"} the correlation between the error terms in the exposure and outcome models.
#'@param progress Logical, indicating whether or not the progress (i.e. the \code{\link{proc.time}} for each \code{Rho}) of the optimization will be output
#'@param ... Additional arguments to be passed on to the \code{maxLik} function. Can be used to set the \code{method} and \code{control} arguments of the \code{maxLik} function.

#'@details
#'The maximization of the log-likelihood is performed using \code{\link[maxLik]{maxLik}}, the default is to use the Newton-Raphson method and an analytic gradient and Hessian.
#'@return \code{coefs.sensmed} returns a list with elements:
#'\item{call}{The matched call}
#' \item{coef}{A matrix with the estimated regression parameters for \code{model.resp} over the range of \code{Rho}. One column per value of \code{Rho}.}
#' \item{sigma.res.resp}{If \code{model.resp} is a linear regression model, the estimated standard deviation of the error term for each \code{Rho}.}
#' \item{sigma.res.expl}{If \code{model.expl} is a linear regression model, the estimated standard deviation of the error term for each \code{Rho}.}
#' \item{Rho}{The sensitivity parameter vector.}
#' \item{expl.coef}{A matrix with the estimated regression parameters for \code{model.expl} over the range of \code{Rho}. One column per value of \code{Rho}.}
#' \item{model.expl}{the original fitted \code{glm} object of \code{model.expl}.}
#' \item{model.resp}{the original fitted \code{glm} object of \code{model.resp}.}
#' \item{X.expl}{The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}}
#' \item{X.resp}{The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}}
#' \item{outc.resp}{The outcome variable of \code{model.resp}.}
#' \item{outc.expl}{The outcome variable of \code{model.expl}.}
#' \item{sigmas}{A list with the estimated covariance matrices for the regression parameters of \code{model.resp} and \code{model.expl} over \code{Rho}.}
#' \item{max.info}{Information about the maximization (whether or not the convergence was successful, \code{message}, \code{method} and number of iterations) for each \code{Rho}, see \code{\link[maxLik]{maxLik}} for more information.}
#' \item{value}{The values of the loglikelihood function for the best set of regression parameters from the optimization for each \code{Rho}, see \code{\link[maxLik]{maxLik}}.}
#' @author Anita Lindmark
#' @references Henningsen, A., Toomet, O. (2011). maxLik: A Package for Maximum Likelihood Estimation in R, \emph{Computational Statistics}, \bold{26(3)}, pp. 443--458.
#' @export
#' @seealso \code{\link{sensmediation}}
#'@examples
#'
#'\dontrun{
#' # Example with data from Riksstroke (the Swedish stroke register)
#'
#' data(RSdata)
#'
#' # Probit mediator and outcome models:
#' m.model <- glm(lowered.consc ~ AF + age.cat + sex, data = RSdata,
#'    family = binomial(link = 'probit'))
#' o.model <- glm(cf.3mo ~ AF + lowered.consc + age.cat + sex, data = RSdata,
#'    family = binomial(link = 'probit'))
#'
#' # Estimation of regression coefficients under different values of Rho
#' # Rho = correlation between error terms in mediator and outcome model:
#' coefs.MY <- coefs.sensmed(model.expl = m.model, model.resp = o.model, Rho = seq(0, 0.5, 0.1))
#' # Outcome model regression coefficients:
#' coefs.MY$coef
#'}
#'
coefs.sensmed <- function(model.expl, model.resp, Rho, progress = TRUE, ...){

  cll <- match.call()

  # Check to make sure 0 is part of Rho
  if(sum(Rho == 0) == 0){
    Rho <- c(0, Rho)
    warning('Note that, 0 was not an element of Rho as recommended. Hence 0 is now added into Rho.')
  }

  Rho <- sort(Rho) # Sort Rho from smallest to largest value
  nrho <- length(Rho)
  i0 <- which(Rho == 0)

  # Selecting which maximum likelihood function to run depending on the types of models input
  if(model.expl$family$link == "probit"){
    if(model.resp$family$link == "probit"){ # Probit model.resp and model.expl
      ML <- ML.bb(model.expl, model.resp, Rho, progress, ...) # Estimates of regression parameters

    }
    if(model.resp$family$family == "gaussian"){ # Linear model.resp, probit model.expl
      ML <- ML.bc(model.expl, model.resp, Rho, progress, ...) # Estimates of regression parameters

    }
  }

  if(model.expl$family$family == "gaussian"){
    if(model.resp$family$link == "probit"){ # Probit model.resp, linear model.expl
      ML <- ML.cb(model.expl, model.resp, Rho, progress, ...) # Estimates of regression parameters
      }
    if(model.resp$family$family == "gaussian"){ # Linear model.resp and model.expl
      ML <- ML.cc(model.expl, model.resp, Rho, progress, ...) # Estimates of regression parameters

    }
  }

  ML$call <- cll

  return(ML)
}

#'Functions for ML estimation of regression parameters for sensitivity analysis
#'
#'Functions for ML estimation of regression parameters for sensitivity analysis for different combinations of exposure, mediator and outcome models. The functions are named according to the convention \code{ML."model.expl type""model.resp type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression. The optimization is performed using
#'\code{\link[maxLik]{maxLik}}. The functions are intended to be called through \code{\link{coefs.sensmed}}, not on their own.
#'
#'@param model.expl Fitted \code{\link{glm}} model object (probit or linear). If sensitivity analysis to mediator-outcome confounding the mediator model. Otherwise the exposure model.
#'@param model.resp Fitted \code{\link{glm}} model object (probit or linear). If sensitivity analysis to exposure-mediator confounding the mediator model. Otherwise the outcome model.
#'@param Rho The sensitivity parameter vector. If \code{type="my"} the correlation between the error terms in the mediator and outcome models. If \code{type="zm"} the correlation between the error terms in the exposure and mediator models. If \code{type="zy"} the correlation between the error terms in the exposure and outcome models.
#'@param progress Logical, indicating whether or not the progress (i.e. the \code{\link{proc.time}} for each \code{Rho}) of the optimization will be output
#'@param ... Additional arguments to be passed on to the \code{maxLik} function. Can be used to set the \code{method} and \code{control} arguments of the \code{maxLik} function.
#'@return A list with elements:
#' \item{coef}{A matrix with the estimated regression parameters for \code{model.resp} over the range of \code{Rho}. One column per value of \code{Rho}.}
#' \item{Rho}{The sensitivity parameter vector.}
#' \item{expl.coef}{A matrix with the estimated regression parameters for \code{model.expl} over the range of \code{Rho}. One column per value of \code{Rho}.}
#' \item{model.expl}{the original fitted \code{glm} object of \code{model.expl}.}
#' \item{model.resp}{the original fitted \code{glm} object of \code{model.resp}.}
#' \item{X.expl}{The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}}
#' \item{X.resp}{The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}}
#' \item{outc.resp}{The outcome variable of \code{model.resp}.}
#' \item{outc.expl}{The outcome variable of \code{model.expl}.}
#' \item{sigma.res.expl}{If \code{model.expl} is linear, a column matrix with the estimated residual standard deviation for \code{model.expl} over the range of \code{Rho}.}
#' \item{sigma.res.resp}{If \code{model.resp} is linear, a column matrix with the estimated residual standard deviation for \code{model.resp} over the range of \code{Rho}.}
#' \item{value}{The values of the -loglikelihood function for the best set of regression parameters from the optimization for each \code{Rho}.}
#' \item{sigmas}{A list with the covariance matrices for the model parameters in \code{model.expl} and \code{model.resp} for each \code{Rho}.}
#' \item{max.info}{Information about the maximization (whether or not the convergence was successful, \code{message}, \code{method} and number of iterations) for each \code{Rho}, see \code{\link[maxLik]{maxLik}} for more information.}
#' @author Anita Lindmark
#' @seealso \code{\link{coefs.sensmed}}, \code{\link[maxLik]{maxLik}}
#' @name ML
NULL


#' @rdname ML
#' @export
ML.bb <- function(model.expl, model.resp, Rho, progress = TRUE, ...){

  X.expl <- as.matrix(stats::model.matrix(model.expl))
  X.resp <- as.matrix(stats::model.matrix(model.resp))

  outc.resp <- model.resp$y
  outc.expl <- as.matrix(model.expl$y)

  nrho <- length(Rho)
  p1 <- length(coef(model.resp))
  p2 <- length(coef(model.expl))
  d <- p1 + p2

  coefs <- matrix(nrow=(p1 + p2), ncol = nrho) # vector with all model parameters
  rownames(coefs) <- c(names(model.resp$coef), names(model.expl$coef))

  dots <- list(...)

  if(is.null(dots$method)) # which optimization method should be used?
    method <- "NR"
  else
    method <- dots$method

  if(is.null(dots$control)) # any control arguments given?
    control <- NULL
  else
    control <- dots$control

  i0 <- which(Rho == 0)
  coefs[, i0] <- c(coef(model.resp), coef(model.expl)) # storing the parameters for Rho = 0

  # Convergence information, matrix and stored values for Rho = 0
  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- model.expl$converged == TRUE & model.resp$converged == TRUE
  max.info$message[i0] <- " "
  max.info$method[i0] <- "glm"
  max.info$iterations[i0] <- NA

  value <- vector(length = nrho) # vector to store values of the log-likelihood
  names(value) <- paste(Rho)
  value[i0] <- LogL.bb(coefs[, i0], Rho = 0, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
                       outc.expl = outc.expl)



  sigma <- list() # List to store the covariance matrices for the model parameters
  sigma[[i0]] <- matrix(0, nrow = d, ncol = d)
  glmVcov <- stats::vcov(model.resp)
  expl.glmVcov <- stats::vcov(model.expl)
  sigma[[i0]][1:p1, 1:p1] <- glmVcov
  sigma[[i0]][(p1 + 1):d, (p1 + 1):d] <- expl.glmVcov
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], dimnames(expl.glmVcov)[[1]])
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  # Optimization for Rho < 0:
  if(sum(Rho < 0) > 0){
    for(i in 1:(i0 - 1)){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.bb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.bb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.bb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i0 - i + 1], method = method, control = control)

      # Storing information from maximization
      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i0 - i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Optimization for Rho > 0:
  if(sum(Rho > 0) > 0){
    for(i in (i0 + 1):nrho){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.bb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.bb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.bb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)


      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i - 1], method = method, control = control)

      # Storing information from maximization
      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i] <- ML$estimate
      value[i]  <-  ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Storing the estimated parameters for model.expl and model.resp over Rho.
  expl.coef  <- as.matrix(coefs[(p1 + 1):(p1 + p2), ])
  if(p2 == 1){
    expl.coef  <- t(expl.coef)
    rownames(expl.coef) <- names(model.expl$coefficients)
  }

  coef <- as.matrix(coefs[1:(p1), ])
  if(p1 == 1){
    coef <- t(coef)
    rownames(coef) <- names(model.resp$coefficients)
  }
  colnames(expl.coef) <- colnames(coef) <- names(sigma) <- paste(Rho)

  max.info <- lapply(max.info, stats::setNames, nm = paste(Rho))

  return(list(coef = coef, Rho = Rho, expl.coef = expl.coef, model.expl = model.expl, model.resp = model.resp, X.expl = X.expl,
              X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, value = value, sigmas = sigma, max.info = max.info))

}

#' @rdname ML
#' @export
ML.bc <- function(model.expl, model.resp, Rho, progress = TRUE, ...){

  X.expl <- as.matrix(stats::model.matrix(model.expl))
  X.resp <- as.matrix(stats::model.matrix(model.resp))

  outc.resp <- model.resp$y
  outc.expl <- as.matrix(model.expl$y)

  nrho <- length(Rho)
  p1 <- length(coef(model.resp))
  p2 <- length(coef(model.expl))

  d <- p1 + p2 + 1

  coefs <- matrix(nrow=(p1 + p2 + 1), ncol = nrho) # vector with all model parameters
  rownames(coefs) <- c(names(model.resp$coef), "sigma.res.resp", names(model.expl$coef))


  dots <- list(...)
  if(is.null(dots$method)) # which optimization method should be used?
    method <- "NR"
  else
    method <- dots$method
  if(is.null(dots$control)) # any control arguments given?
    control <- NULL
  else
    control <- dots$control

  i0 <- which(Rho == 0)
  coefs[, i0] <- c(coef(model.resp), sqrt(summary(model.resp)$dispersion), coef(model.expl)) # storing the parameters for Rho = 0

  # Convergence information, matrix and stored values for Rho = 0
  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- model.expl$converged == TRUE & model.resp$converged == TRUE
  max.info$message[i0] <- " "
  max.info$method[i0] <- "glm"
  max.info$iterations[i0] <- NA


  value <- vector(length = nrho) # vector to store values of the log-likelihood
  names(value) <- paste(Rho)
  value[i0] <- LogL.bc(coefs[, i0], Rho = 0, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
                       outc.expl = outc.expl)

  sigma <- list() # List to store the covariance matrices for the model parameters
  sigma[[i0]] <- matrix(0, nrow = d, ncol = d)
  glmVcov <- stats::vcov(model.resp)
  expl.glmVcov <- stats::vcov(model.expl)
  sigma[[i0]][1:p1, 1:p1] <- glmVcov
  sigma[[i0]][(p1 + 2):d, (p1 + 2):d] <- expl.glmVcov
  sigma[[i0]][(p1 + 1), (p1 + 1)] <- summary(model.resp)$dispersion^2/(2*model.resp$df.residual)
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], "sigma.res.resp", dimnames(expl.glmVcov)[[1]])
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  # Optimization for Rho < 0:
  if(sum(Rho < 0) > 0){
    for(i in 1:(i0 - 1)){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.bc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.bc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.bc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i0 - i + 1], method = method)

      # Storing information from maximization
      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i0 - i]] <- -solve(ML$hessian)


      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Optimization for Rho > 0:
  if(sum(Rho > 0) > 0){
    for(i in (i0 + 1):nrho){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.bc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.bc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.bc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i - 1], method = method)

      # Storing information from maximization
      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i] <- ML$estimate
      value[i]  <-  ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Storing the estimated parameters for model.expl and model.resp over Rho.
  expl.coef <- as.matrix(coefs[(p1 + 2):(p1 + p2 + 1), ])
  if(p2 == 1){
    expl.coef  <- t(expl.coef)
    rownames(expl.coef) <- names(model.expl$coefficients)
  }
  sigma.res.resp <- as.matrix(coefs[p1 + 1, ])
  coef <- as.matrix(coefs[1:p1, ])
  if(p1 == 1){
    coef <- t(coef)
    rownames(coef) <- names(model.resp$coefficients)
  }
  colnames(expl.coef) <- rownames(sigma.res.resp) <- colnames(coef) <- names(sigma) <- paste(Rho)

  max.info <- lapply(max.info, stats::setNames, nm = paste(Rho))

  return(list(coef = coef, Rho = Rho, expl.coef = expl.coef, model.expl = model.expl,
              model.resp = model.resp, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl,
              value = value, sigma.res.resp = sigma.res.resp, sigmas = sigma, max.info = max.info))


}

#' @rdname ML
#' @export
ML.cb <- function(model.expl, model.resp, Rho, progress = TRUE, ...){

  X.expl <- as.matrix(stats::model.matrix(model.expl))
  X.resp <- as.matrix(stats::model.matrix(model.resp))

  outc.resp <- model.resp$y
  outc.expl <- as.matrix(model.expl$y)

  nrho <- length(Rho)
  p1 <- length(coef(model.resp))
  p2 <- length(coef(model.expl)) + 1

  d <- p1 + p2 - 1

  coefs <- matrix(nrow = (p1 + p2), ncol = nrho) # vector with all model parameters
  rownames(coefs) <- c(names(model.resp$coef), names(model.expl$coef), "sigma.res.expl")

  dots <- list(...)
  if(is.null(dots$method)) # which optimization method should be used?
    method <- "NR"
  else
    method <- dots$method
  if(is.null(dots$control)) # any control arguments given?
    control <- NULL
  else
    control <- dots$control

  i0 <- which(Rho == 0)
  coefs[, i0] <- c(coef(model.resp), coef(model.expl), sqrt(summary(model.expl)$dispersion)) # storing the parameters for Rho = 0

  # Convergence information, matrix and stored values for Rho = 0
  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- model.expl$converged == TRUE & model.resp$converged == TRUE
  max.info$message[i0] <- " "
  max.info$method[i0] <- "glm"
  max.info$iterations[i0] <- NA

  value <- vector(length = nrho) # vector to store values of the log-likelihood
  names(value) <- paste(Rho)
  value[i0] <- LogL.cb(coefs[, i0], Rho = 0, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
                       outc.expl = outc.expl)

  sigma <- list() # List to store the covariance matrices for the model parameters
  sigma[[i0]] <- matrix(0, nrow = d + 1, ncol = d + 1)
  glmVcov <- stats::vcov(model.resp)
  expl.glmVcov <- stats::vcov(model.expl)
  sigma[[i0]][1:p1, 1:p1] <- glmVcov
  sigma[[i0]][(p1 + 1):d, (p1 + 1):d] <- expl.glmVcov
  sigma[[i0]][(d + 1), (d + 1)] <- summary(model.expl)$dispersion/(2*model.expl$df.residual)
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], dimnames(expl.glmVcov)[[1]], "sigma.res.expl")
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  # Optimization for Rho < 0:
  if(sum(Rho < 0) > 0){
    for(i in 1:(i0 - 1)){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.cb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.cb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.cb(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i0 - i + 1], method = method)

      # Storing information from maximization
      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i0 - i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Optimization for Rho > 0:
  if(sum(Rho > 0) > 0){
    for(i in (i0 + 1):nrho){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.cb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.cb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.cb(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)


      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i - 1], method = method)

      # Storing information from maximization
      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i] <- ML$estimate
      value[i]  <-  ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Storing the estimated parameters for model.expl and model.resp over Rho.
  expl.coef <- as.matrix(coefs[(p1 + 1):(p1 + p2 - 1), ])
  if(p2 == 2){
    expl.coef  <- t(expl.coef)
    rownames(expl.coef) <- names(model.expl$coefficients)
  }
  sigma.res.expl <- as.matrix(coefs[p1 + p2, ])
  coef <- as.matrix(coefs[1:(p1), ])
  if(p1 == 1){
    coef <- t(coef)
    rownames(coef) <- names(model.resp$coefficients)
  }
  colnames(expl.coef) <- rownames(sigma.res.expl) <- colnames(coef) <- names(sigma) <- paste(Rho)

  max.info <- lapply(max.info, stats::setNames, nm = paste(Rho))

  return(list(coef = coef, Rho = Rho, expl.coef = expl.coef, model.expl = model.expl,
              model.resp = model.resp, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl,
              value = value, sigma.res.expl = sigma.res.expl, sigmas = sigma, max.info = max.info))

}


#' @rdname ML
#' @export
ML.cc <- function(model.expl, model.resp, Rho, progress = TRUE, ...){

  X.expl <- as.matrix(stats::model.matrix(model.expl))
  X.resp <- as.matrix(stats::model.matrix(model.resp))

  outc.resp <- model.resp$y
  outc.expl <- as.matrix(model.expl$y)

  nrho <- length(Rho)
  p1 <- length(coef(model.resp)) + 1
  p2 <- length(coef(model.expl)) + 1

  coefs <- matrix(nrow = (p1  +p2), ncol = nrho) # vector with all model parameters
  rownames(coefs) <- c(names(model.resp$coef), "sigma.res.resp", names(model.expl$coef), "sigma.res.expl")

  dots <- list(...)
  if(is.null(dots$method)) # which optimization method should be used?
    method <- "NR"
  else
    method <- dots$method
  if(is.null(dots$control)) # any control arguments given?
    control <- NULL
  else
    control <- dots$control

  i0 <- which(Rho == 0)
  coefs[, i0] <- c(coef(model.resp), sqrt(summary(model.resp)$dispersion), coef(model.expl),
                   sqrt(summary(model.expl)$dispersion)) # storing the parameters for Rho = 0

  # Convergence information, matrix and stored values for Rho = 0
  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- model.expl$converged == TRUE & model.resp$converged == TRUE
  max.info$message[i0] <- " "
  max.info$method[i0] <- "glm"
  max.info$iterations[i0] <- NA

  value <- vector(length = nrho) # vector to store values of the log-likelihood
  names(value) <- paste(Rho)
  value[i0] <- LogL.cc(coefs[, i0], Rho = 0, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
                       outc.expl = outc.expl)

  sigma <- list() # List to store the covariance matrices for the model parameters
  sigma[[i0]] <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)
  glmVcov <- stats::vcov(model.resp)
  expl.glmVcov <- stats::vcov(model.expl)
  sigma[[i0]][1:(p1 - 1), 1:(p1 - 1)] <- glmVcov
  sigma[[i0]][(p1 + 1):(p1 + p2 - 1), (p1 + 1):(p1 + p2 - 1)] <- expl.glmVcov
  sigma[[i0]][(p1), (p1)] <- summary(model.resp)$dispersion/(2*model.resp$df.residual)
  sigma[[i0]][(p1 + p2), (p1 + p2)] <- summary(model.expl)$dispersion/(2*model.expl$df.residual)
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], "sigma.res.resp", dimnames(expl.glmVcov)[[1]], "sigma.res.expl")
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  # Optimization for Rho < 0:
  if(sum(Rho < 0) > 0){
    for(i in 1:(i0-1)){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.cc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.cc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.cc(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i0 - i + 1], method = method)

      # Storing information from maximization
      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i0 - i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Optimization for Rho > 0:
  if(sum(Rho > 0) > 0){
    for(i in (i0 + 1):nrho){

      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }

      # Functions for use by maxLik, log-likelihood and analytic gradient and hessian
      f <- function(par)
        LogL.cc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      g <- function(par)
        grr.cc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      h <- function(par)
        hess.cc(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl)

      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, hess = h, start = coefs[, i - 1], method = method)

      # Storing information from maximization
      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations

      # Storing parameters and log-likelihood value from maximization
      coefs[, i] <- ML$estimate
      value[i]  <-  ML$maximum

      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  # Storing the estimated parameters for model.expl and model.resp over Rho.
  expl.coef <-as.matrix(coefs[(p1 + 1):(p1 + p2 - 1), ])
  if(p2 == 2){
    expl.coef  <- t(expl.coef)
    rownames(expl.coef) <- names(model.expl$coefficients)
  }
  sigma.res.expl <- as.matrix(coefs[p1 + p2, ])
  coef <- as.matrix(coefs[1:(p1 - 1), ])
  if(p1 == 2){
    coef <- t(coef)
    rownames(coef) <- names(model.resp$coefficients)
  }
  sigma.res.resp <- as.matrix(coefs[p1, ])
  colnames(expl.coef) <- rownames(sigma.res.expl) <- rownames(sigma.res.resp) <- colnames(coef) <- names(sigma) <- paste(Rho)

  max.info <- lapply(max.info, stats::setNames, nm = paste(Rho))

  return(list(coef = coef, Rho = Rho, expl.coef = expl.coef, model.expl = model.expl,
              model.resp = model.resp, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
              outc.expl = outc.expl, value = value, sigma.res.resp = sigma.res.resp, sigma.res.expl = sigma.res.expl,
              sigmas = sigma, max.info = max.info))

}

#' Implementation of loglikelihood functions for ML estimation of regression parameters
#'
#'Implementation of loglikelihood functions for ML estimation of regression parameters for different combinations of
#'exposure, mediator and outcome models. The functions are named according to the convention \code{LogL."model.expl type""model.resp type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression.

#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}
#'@param outc.resp The outcome of \code{model.resp}, a vector.
#'@param outc.expl The outcome of \code{model.expl}, a column matrix.
#' @seealso \code{\link{coefs.sensmed}}, \code{\link[maxLik]{maxLik}}
#' @name LogL
NULL


#' @rdname LogL
#' @export
LogL.bb <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  coef.expl <- as.matrix(par[(d.resp + 1):(length(par))])

  # Components to be used in the log-likelihood
  w.resp <- (2*outc.resp - 1)*X.resp%*%coef.resp
  w.expl <- X.expl%*%coef.expl
  Rhos <- (2*outc.resp - 1)*Rho
  n <- length(outc.resp)

  # Calculating the terms of the log-likelihood:
  terms <- numeric(0)
  for(i in 1:n){
    if(outc.expl[i] == 1)
      terms[i] <- log(mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp[i],w.expl[i]), mean = c(0, 0),
                              corr = rbind(c(1, Rhos[i]), c(Rhos[i], 1))))
    else
      terms[i] <- log(mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp[i], -w.expl[i]), mean = c(0, 0),
                              corr = rbind(c(1, -Rhos[i]), c(-Rhos[i], 1))))
  }


  logl <- sum(terms)

  return(logl)
}



#' @rdname LogL
#' @export
LogL.bc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par))])
  sigma.res.resp <- par[(d.resp + 1)]

  # Components to be used in the log-likelihood
  n <- length(outc.resp)
  w1 <- X.expl%*%expl.coef
  w2 <- (outc.resp - X.resp%*%coef.resp)/sigma.res.resp

  logl <- -n*log(sigma.res.resp) +
    sum(log(stats::pnorm((2*outc.expl - 1)*(w1 + Rho*w2)/sqrt(1 - Rho^2))) +
          log(stats::dnorm(w2)))

  return(logl)
}


#' @rdname LogL
#' @export
LogL.cb <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 1):(length(par) - 1)])
  sigma.res.expl <- par[length(par)]

  # Components to be used in the log-likelihood
  n <- length(outc.resp)
  w1 <- X.resp%*%coef.resp
  w2 <- (outc.expl - X.expl%*%expl.coef)/sigma.res.expl

  logl <- -n*log(sigma.res.expl) + sum(log(stats::pnorm((2*outc.resp - 1)*(w1 + Rho*w2)/sqrt(1 - Rho^2))) +
                                         log(stats::dnorm(w2)))

  return(logl)
}


#' @rdname LogL
#' @export
LogL.cc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp+2):(length(par)-1)] )
  sigma.res.resp <- par[d.resp+1]
  sigma.res.expl <- par[length(par)]

  # Components to be used in the log-likelihood
  Q1 <- outc.expl - X.expl%*%expl.coef
  Q2 <- outc.resp - X.resp%*%coef.resp
  Qs <- cbind(Q1, Q2)
  Sigma <- cbind(c(sigma.res.expl^2, Rho*sigma.res.expl*sigma.res.resp),
                 c(Rho*sigma.res.expl*sigma.res.resp, sigma.res.resp^2))

  logl <- sum(log(mvtnorm::dmvnorm(Qs, sigma = Sigma)))

  return(logl)
}


#' Analytic gradients of the loglikelihood functions for ML estimation of regression parameters
#'
#'Implementation of the analytic gradients of the loglikelihood functions for ML estimation of regression parameters for different combinations of
#'exposure, mediator and outcome models. The functions are named according to the convention \code{grr."model.expl type""model.resp type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression.

#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}
#'@param outc.resp The outcome of \code{model.resp}, a vector.
#'@param outc.expl The outcome of \code{model.expl}, a column matrix.
#' @seealso \code{\link{coefs.sensmed}}, \code{\link[maxLik]{maxLik}}
#' @name grr
NULL


#' @rdname grr
#' @export
grr.bb <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- par[1:d.resp]
  coef.expl <- par[(d.resp + 1):(length(par))]

  # Components to be used in the calculation of the gradients
  q <- 2*outc.resp - 1
  w.resp1 <- tcrossprod(coef.resp, X.resp)
  w.resp2 <- q*w.resp1
  w.expl <- tcrossprod(coef.expl, X.expl)
  Rhos <- q*Rho
  n <- length(outc.resp)

  Phi2 <- numeric(0)
  for(i in 1:n){
    if(outc.expl[i] == 1)
      Phi2[i] <- mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp2[i], w.expl[i]), mean = c(0, 0), corr = rbind(c(1, Rhos[i]), c(Rhos[i], 1)))
    else
      Phi2[i] <- mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp2[i], -w.expl[i]), mean = c(0, 0), corr = rbind(c(1, -Rhos[i]), c(-Rhos[i], 1)))
  }
  #

  # Calculating the gradients of the log-likelihood wrt to the parameters of model.resp and model.expl
  gr.d <- ifelse(outc.expl == 1, (q*stats::dnorm(w.resp1)*stats::pnorm((w.expl - Rho*(w.resp1))/sqrt(1 - Rho^2))/Phi2),
                 (q*stats::dnorm(w.resp1)*stats::pnorm((-w.expl + Rho*(w.resp1))/sqrt(1 - Rho^2))/Phi2))

  gr.coef.resp <-  crossprod(gr.d, X.resp) # Gradient of the log-likelihood wrt the parameters of model.resp

  gr.b <- ifelse(outc.expl == 1, (stats::dnorm(w.expl)*stats::pnorm(q*(w.resp1 - Rho*(w.expl))/sqrt(1 - Rho^2))/Phi2),
                 -(stats::dnorm(w.expl)*stats::pnorm(q*(w.resp1 - Rho*(w.expl))/sqrt(1 - Rho^2))/Phi2))

  gr.coef.expl <-  crossprod(gr.b, X.expl) # Gradient of the log-likelihood wrt the parameters of model.expl

  return(c(gr.coef.resp, gr.coef.expl))
}

#' @rdname grr
#' @export
grr.bc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par))])
  sigma.res.resp <- par[(d.resp + 1)]

  # Components to be used in the calculation of the gradients
  n <- length(outc.resp)
  w1 <- X.expl%*%expl.coef
  w2 <- (outc.resp - X.resp%*%coef.resp)/sigma.res.resp
  q <- 2*outc.expl - 1
  A <- q*(w1 + Rho*w2)/sqrt(1 - Rho^2)
  ratioA <- stats::dnorm(A)/stats::pnorm(A)

  # Calculating the gradients of the log-likelihood wrt to the parameters of model.resp and model.expl
  gr.expl.coef <- crossprod(q/sqrt(1 - Rho^2)*ratioA, X.expl) # Gradient of the log-likelihood wrt the parameters of model.expl

  gr.coef.resp <- crossprod(-q*Rho/(sigma.res.resp*sqrt(1-Rho^2))*ratioA + w2/sigma.res.resp, X.resp) # Gradient of the log-likelihood wrt the parameters of model.resp

  gr.sigma <- -n/sigma.res.resp  +  sum(w2^2/sigma.res.resp  -
                                          q*Rho*as.vector(w2/(sigma.res.resp*sqrt(1 - Rho^2)))*as.vector(ratioA)) # Gradient of the log-likelihood wrt the error term standard deviation of model.resp

  return(c(gr.coef.resp, gr.sigma, gr.expl.coef))

}

#' @rdname grr
#' @export
grr.cb <-  function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 1):(length(par) - 1)])
  sigma.res.expl <- par[length(par)]

  # Components to be used in the calculation of the gradients
  n <- length(outc.resp)
  w1 <-  X.resp%*%coef.resp
  w2 <- (outc.expl - X.expl%*%expl.coef)/sigma.res.expl
  q <- 2*outc.resp - 1
  B <- q*(w1 + Rho*w2)/sqrt(1 - Rho^2)
  ratioB <- stats::dnorm(B)/stats::pnorm(B)

  # Calculating the gradients of the log-likelihood wrt to the parameters of model.resp and model.expl
  gr.coef.resp <- crossprod(q/sqrt(1 - Rho^2)*ratioB, X.resp) # Gradient of the log-likelihood wrt the parameters of model.resp

  gr.expl.coef <- crossprod(-q*Rho/(sigma.res.expl*sqrt(1 - Rho^2))*ratioB + w2/sigma.res.expl, X.expl) # Gradient of the log-likelihood wrt the parameters of model.expl

  gr.sigma <- -n/sigma.res.expl + sum(w2^2/sigma.res.expl -
                                        q*Rho*as.vector(w2/(sigma.res.expl*sqrt(1 - Rho^2)))*as.vector(ratioB)) # Gradient of the log-likelihood wrt the error term standard deviation of model.expl

  return(c(gr.coef.resp, gr.expl.coef, gr.sigma))
}

#' @rdname grr
#' @export
grr.cc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)])
  sigma.res.resp <- par[d.resp + 1]
  sigma.res.expl <- par[length(par)]

  # Components to be used in the calculation of the gradients
  n <- length(outc.resp)
  q1 <- outc.expl - X.expl%*%expl.coef
  q2 <- outc.resp - X.resp%*%coef.resp

  # Calculating the gradients of the log-likelihood wrt to the parameters of model.resp and model.expl
  gr.expl.coef <- colSums(as.vector(q1/sigma.res.expl -
                                      Rho*q2/sigma.res.resp)*X.expl)*1/((1 - Rho^2)*sigma.res.expl) # Gradient of the log-likelihood wrt the parameters of model.expl

  gr.coef.resp <- colSums(as.vector(q2/sigma.res.resp -
                                      Rho*q1/sigma.res.expl)*X.resp)/((1 - Rho^2)*sigma.res.resp) # Gradient of the log-likelihood wrt the parameters of model.resp

  gr.sigma.b <- -n/sigma.res.expl + 1/(sigma.res.expl^2*(1 - Rho^2))*sum(q1^2/sigma.res.expl -
                                                                           Rho*q2*q1/sigma.res.resp) # Gradient of the log-likelihood wrt the error term standard deviation of model.expl

  gr.sigma.t <- -n/sigma.res.resp + 1/(sigma.res.resp^2*(1 - Rho^2))*sum(q2^2/sigma.res.resp -
                                                                           Rho*q2*q1/sigma.res.expl) # Gradient of the log-likelihood wrt the error term standard deviation of model.resp

  return(c(gr.coef.resp, gr.sigma.t, gr.expl.coef, gr.sigma.b))

}

#' Analytic Hessians of the loglikelihood functions for ML estimation of regression parameters
#'
#'Implementation of the analytic Hessians of the loglikelihood functions for ML estimation of regression parameters for different combinations of
#'exposure, mediator and outcome models. The functions are named according to the convention \code{hess."model.expl type""model.resp type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression.

#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}
#'@param outc.resp The outcome of \code{model.resp}, a vector.
#'@param outc.expl The outcome of \code{model.expl}, a column matrix.
#' @seealso \code{\link{coefs.sensmed}}, \code{\link[maxLik]{maxLik}}
#' @name hess
NULL


#' @rdname hess
#' @export
hess.bb <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- par[1:d.resp]
  coef.expl <- par[(d.resp + 1):(length(par))]

  # Components to be used in the calculation of the hessians
  w.expl <- tcrossprod(coef.expl, X.expl)
  w.resp1 <- tcrossprod(coef.resp, X.resp)
  q <- 2*outc.resp - 1
  w.resp2 <- q*w.resp1
  Rhos <- q*Rho
  n <- length(outc.resp)

  Phi2 <- numeric(0)
  for(i in 1:n){
    if(outc.expl[i] == 1)
      Phi2[i] <- mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp2[i], w.expl[i]), mean = c(0, 0), corr = rbind(c(1, Rhos[i]), c(Rhos[i], 1)))
    else
      Phi2[i] <- mvtnorm::pmvnorm(lower = -Inf, upper = c(w.resp2[i], -w.expl[i]), mean = c(0, 0), corr = rbind(c(1, -Rhos[i]), c(-Rhos[i], 1)))
  }

  pnorm1 <- stats::pnorm((w.resp2 - q*Rho*(w.expl))/sqrt(1 - Rho^2))
  dnorm1 <- stats::dnorm((w.resp2 - q*Rho*(w.expl))/sqrt(1 - Rho^2))
  dnorm.expl <- stats::dnorm(w.expl)
  dnorm.resp <- stats::dnorm(w.resp1)
  pnorm2 <- stats::pnorm((w.expl - Rho*(w.resp1))/sqrt(1 - Rho^2))
  dnorm2 <- stats::dnorm((w.expl - Rho*(w.resp1))/sqrt(1 - Rho^2))
  pnorm3 <- stats::pnorm((-w.expl + Rho*(w.resp1))/sqrt(1 - Rho^2))
  dnorm3 <- stats::dnorm((-w.expl + Rho*(w.resp1))/sqrt(1 - Rho^2))
  ##

  # Calculating the hessians of the log-likelihood wrt to the parameters of model.resp and model.expl
  secder.expl <- ifelse(outc.expl == 1,
                        -dnorm.expl/Phi2*(pnorm1^2/Phi2*dnorm.expl + w.expl*pnorm1 + dnorm1*q*Rho/sqrt(1 - Rho^2)),
                        -dnorm.expl/Phi2*(pnorm1^2/Phi2*dnorm.expl - w.expl*pnorm1 - dnorm1*q*Rho/sqrt(1 - Rho^2)))

  hess.expl <- crossprod(c(secder.expl)*X.expl, X.expl) # Second partial derivative wrt the parameters of model.expl

  secder.resp <- ifelse(outc.expl == 1,
                        -q*dnorm.resp/Phi2*(q*pnorm2^2/Phi2*dnorm.resp + w.resp1*pnorm2 + dnorm2*Rho/sqrt(1 - Rho^2)),
                        q*dnorm.resp/Phi2*(-q*pnorm3^2/Phi2*dnorm.resp - w.resp1*pnorm3 + dnorm3*Rho/sqrt(1 - Rho^2)))

  hess.resp <- crossprod(c(secder.resp)*X.resp, X.resp) # Second partial derivative wrt the parameters of model.resp

  secder.er <- ifelse(outc.expl==1,
                      q*dnorm.expl/Phi2*(-dnorm.resp/Phi2*pnorm1*pnorm2+dnorm1/sqrt(1-Rho^2)),
                      -q*dnorm.expl/Phi2*(-dnorm.resp/Phi2*pnorm1*pnorm3+dnorm1/sqrt(1-Rho^2)))

  hess.er <- crossprod(c(secder.er)*X.expl, X.resp) # Mixed partial derivative wrt the parameters of model.expl and model.resp


  hess <- rbind(cbind(hess.resp, t(hess.er)), cbind(hess.er, hess.expl))

  return(hess)
}

#' @rdname hess
#' @export
hess.bc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par))])
  sigma.res.resp <- par[(d.resp + 1)]

  # Components to be used in the calculation of the hessians
  n <- length(outc.resp)
  w1 <- X.expl%*%expl.coef
  w2 <- (outc.resp - X.resp%*%coef.resp)/sigma.res.resp
  q <- 2*outc.expl - 1
  A <- q*(w1 + Rho*w2)/sqrt(1 - Rho^2)
  ratioA <- stats::dnorm(A)/stats::pnorm(A)

  # Calculating the hessians of the log-likelihood wrt to the parameters of model.resp and model.expl
  secder.expl <- -ratioA/(1 - Rho^2)*(ratioA + A)

  hess.expl <- crossprod(c(secder.expl)*X.expl, X.expl) # Second partial derivative wrt the parameters of model.expl

  secder.er <- ratioA*Rho/(sigma.res.resp*(1 - Rho^2))*(ratioA + A)

  hess.er <- crossprod(c(secder.er)*X.expl, X.resp) # Mixed partial derivative wrt the parameters of model.expl and model.resp

  secder.resp <- -1/sigma.res.resp^2*(ratioA*Rho^2/(1 - Rho^2)*(ratioA + A) + 1)

  hess.resp <- crossprod(c(secder.resp)*X.resp, X.resp) # Second partial derivative wrt the parameters of model.resp

  secder.es <- Rho*w2/(sigma.res.resp*(1 - Rho^2))*ratioA*(ratioA + A)

  hess.es <- crossprod(c(secder.es), X.expl) # Mixed partial derivative wrt the parameters of model.expl and the error term standard deviation of model.resp

  secder.sr <- 1/(sigma.res.resp^2)*(q*Rho/(sqrt(1 - Rho^2))*ratioA -
                                       Rho^2*w2/(1 - Rho^2)*(ratioA^2 + ratioA*A)  - 2*w2)

  hess.sr <- crossprod(c(secder.sr), X.resp) # Mixed partial derivative wrt the parameters of model.resp and the error term standard deviation of model.resp

  secder.ss <- ratioA/(sigma.res.resp^2)*(-ratioA*Rho^2*w2^2/(1 - Rho^2) -
                                            A*Rho^2*w2^2/(1 - Rho^2) + 2*q*Rho*w2/sqrt(1 - Rho^2)) -
    3*w2^2/(sigma.res.resp^2)

  hess.ss <- sum(secder.ss) + n/(sigma.res.resp^2) # Second partial derivative wrt the error term standard deviation of model.resp

  hess <- rbind(cbind(hess.resp, t(hess.sr), t(hess.er)), cbind(hess.sr, hess.ss, hess.es),
                cbind(hess.er, t(hess.es), hess.expl))

}

#' @rdname hess
#' @export
hess.cb <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 1):(length(par) - 1)])
  sigma.res.expl <- par[length(par)]

  # Components to be used in the calculation of the hessians
  n <- length(outc.resp)
  w1 <- X.resp%*%coef.resp
  w2 <- (outc.expl - X.expl%*%expl.coef)/sigma.res.expl
  q <- 2*outc.resp - 1
  B <- q*(w1 + Rho*w2)/sqrt(1 - Rho^2)
  ratioB <- stats::dnorm(B)/stats::pnorm(B)


  # Calculating the hessians of the log-likelihood wrt to the parameters of model.resp and model.expl
  secder.expl <- -1/(sigma.res.expl^2)*(ratioB*Rho^2/(1 - Rho^2)*(ratioB + B) + 1)

  hess.expl <- crossprod(c(secder.expl)*X.expl, X.expl) # Second partial derivative wrt the parameters of model.expl

  secder.er <- ratioB*Rho/(sigma.res.expl*(1 - Rho^2))*(ratioB + B)

  hess.er <- crossprod(c(secder.er)*X.expl, X.resp) # Mixed partial derivative wrt the parameters of model.expl and model.resp

  secder.resp <- -ratioB/(1 - Rho^2)*(ratioB + B)

  hess.resp <- crossprod(c(secder.resp)*X.resp, X.resp) # Second partial derivative wrt the parameters of model.resp

  secder.sr <- Rho*w2/(sigma.res.expl*(1 - Rho^2))*ratioB*(ratioB + B)

  hess.sr <- crossprod(c(secder.sr), X.resp) # Mixed partial derivative wrt the parameters of model.resp and the error term standard deviation of model.expl

  secder.es <- 1/(sigma.res.expl^2)*(q*Rho/(sqrt(1 - Rho^2))*ratioB -
                                       Rho^2*w2/(1 - Rho^2)*(ratioB^2 + ratioB*B)  - 2*w2)

  hess.es <- crossprod(c(secder.es), X.expl) # Mixed partial derivative wrt the parameters of model.expl and the error term standard deviation of model.expl

  secder.ss <- ratioB/(sigma.res.expl^2)*(-ratioB*Rho^2*w2^2/(1 - Rho^2) -
                                            B*Rho^2*w2^2/(1 - Rho^2) + 2*q*Rho*w2/sqrt(1 - Rho^2)) -
    3*w2^2/(sigma.res.expl^2)

  hess.ss <- sum(secder.ss) + n/(sigma.res.expl^2) # Second partial derivative wrt the error term standard deviation of model.expl

  hess <- rbind(cbind(hess.resp, t(hess.er), t(hess.sr)), cbind(hess.er, hess.expl, t(hess.es)),
                cbind(hess.sr, hess.es, hess.ss))

}

#' @rdname hess
#' @export
hess.cc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl){

  # Separating the coefficients from model.resp and model.expl
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)])
  sigma.res.resp <- par[d.resp + 1]
  sigma.res.expl <- par[length(par)]

  # Components to be used in the calculation of the hessians
  n <- length(outc.resp)
  q1 <- outc.expl - X.expl%*%expl.coef
  q2 <- outc.resp - X.resp%*%coef.resp

  # Calculating the hessians of the log-likelihood wrt to the parameters of model.resp and model.expl
  hess.expl <- -crossprod(1/((1 - Rho^2)*sigma.res.expl^2)*X.expl, X.expl) # Second partial derivative wrt the parameters of model.expl

  hess.er <- crossprod(Rho/((1 - Rho^2)*sigma.res.expl*sigma.res.resp)*X.expl, X.resp) # Mixed partial derivative wrt the parameters of model.expl and model.resp

  hess.resp <- -crossprod(1/((1 - Rho^2)*sigma.res.resp^2)*X.resp, X.resp) # Second partial derivative wrt the parameters of model.resp

  secder.see <- 1/((1 - Rho^2)*sigma.res.expl^2)*(-2*q1/sigma.res.expl + Rho*q2/sigma.res.resp)

  hess.see <- crossprod(c(secder.see), X.expl) # Mixed partial derivative wrt the parameters of model.expl and the error term standard deviation of model.expl

  secder.srr <- 1/((1 - Rho^2)*sigma.res.resp^2)*(-2*q2/sigma.res.resp + Rho*q1/sigma.res.expl)

  hess.srr <- crossprod(c(secder.srr), X.resp) # Mixed partial derivative wrt the parameters of model.resp and the error term standard deviation of model.resp

  secder.ser <- Rho/((1 - Rho^2)*sigma.res.expl^2*sigma.res.resp)*q1

  hess.ser <- crossprod(c(secder.ser), X.resp) # Mixed partial derivative wrt the parameters of model.resp and the error term standard deviation of model.expl

  secder.sre <- Rho/((1 - Rho^2)*sigma.res.expl*sigma.res.resp^2)*q2

  hess.sre <- crossprod(c(secder.sre), X.expl) # Mixed partial derivative wrt the parameters of model.expl and the error term standard deviation of model.resp

  secder.sese <- -2/((1 - Rho^2)*sigma.res.expl^3)*(q1^2/sigma.res.expl - Rho*q1*q2/sigma.res.resp) -
    1/((1 - Rho^2)*sigma.res.expl^2)*(q1^2/(sigma.res.expl^2)   )

  hess.sese <- sum(secder.sese) + n/(sigma.res.expl^2) # Second partial derivative wrt the error term standard deviation of model.expl

  secder.srsr <- -2/((1 - Rho^2)*sigma.res.resp^3)*(q2^2/sigma.res.resp - Rho*q1*q2/sigma.res.expl) -
    1/((1 - Rho^2)*sigma.res.resp^2)*(q2^2/(sigma.res.resp^2) )

  hess.srsr <- sum(secder.srsr) + n/(sigma.res.resp^2) # Second partial derivative wrt the error term standard deviation of model.resp

  hess.sesr <- 1/((1 - Rho^2)*sigma.res.expl^2)*sum(Rho*q1*q2/(sigma.res.resp^2)) # Mixed partial derivative wrt the error term standard deviation of model.expl and the error term standard deviation of model.resp

  hess <- rbind(cbind(hess.resp, t(hess.srr), t(hess.er),t(hess.ser)),
                cbind(hess.srr, hess.srsr, hess.sre, hess.sesr),
                cbind(hess.er, t(hess.sre), hess.expl, t(hess.see)),
                cbind(hess.ser, hess.sesr, hess.see, hess.sese))
}




