#'@include coefs.sensmed.R calc.effects.R
NULL
#'Estimate natural direct and indirect effects based on parametric regression models and perform sensitivity analysis
#'
#'Function to estimate the natural direct and indirect effects based on parametric regression models. Standard errors for the effects are calculated using the delta method.
#'The function also gives sensitivity analysis results for unobserved confounding. Implements methods introduced in Lindmark, de Luna and Eriksson (2018).
#'
#'@param med.model Fitted \code{\link{glm}} model object representing the mediator model at the basis of the estimation (see Details for more information).
#'@param out.model Fitted \code{\link{glm}} model object representing the outcome model at the basis of the estimation (see Details for more information).
#'@param type the type of confounding for which the sensitivity analysis is to be performed. \code{type="my"}, the default, corresponds to unobserved mediator-outcome confounding, \code{type="zm"} to exposure-mediator confounding and \code{type="zy"} to exposure-outcome confounding.
#'@param exp.model Fitted \code{\link{glm}} model object representing the exposure model. Should be provided if \code{type="zm"} or \code{type="zy"}.
#'@param exp.name A character string indicating the name of the exposure variable used in the models. Needs to match the name of the exposure found in the output from the fitted glm-models (this is especially important to check for exposures of class \code{factor}).
#'@param med.name A character string indicating the name of the mediator used in the models. Needs to match the name of the mediator found in the output from the outcome glm-model (this is especially important to check for mediators of class \code{factor}).
#'@param Rho The sensitivity parameter vector. If \code{type="my"} the correlation between the error terms in the mediator and outcome models. If \code{type="zm"} the correlation between the error terms in the exposure and mediator models. If \code{type="zy"} the correlation between the error terms in the exposure and outcome models.
#'@param progress Logical, indicating whether or not the progress (i.e. the \code{\link{proc.time}} for each \code{Rho}) of the optimization will be output
#'@param conf.level the confidence level to be used for confidence intervals and uncertainty intervals.
#'@param covariates if conditional effects are to be estimated the named list of covariate values (see Details). Covariates not specified are marginalized over.
#'@param alt.decomposition logical indicating whether or not alternative definitions of the direct and indirect effects should be used (see Details).
#'@param exp.value value of the exposure variable used as the exposure condition, default is 1.
#'@param control.value value of the exposure variable used as the control (unexposed) condition, default is 0.
#'@param covariance,med.full,out.full,all.interactions arguments used in previous versions of the package that are now deprecated.
#'@param ... Additional arguments to be passed on to the \code{maxLik} function. Can be used to set the \code{method} and \code{control} arguments of the \code{maxLik} function (see \code{\link{coefs.sensmed}}).

#'@details
#'To obtain the ML estimates of the regression parameters used to calculate mediation effects and perform sensitivity analysis
#'\code{sensmediation} calls \code{\link{coefs.sensmed}}. The maximization of the log-likelihood is performed using
#'\code{\link[maxLik]{maxLik}}, the default is to use the Newton-Raphson method and an analytic gradient and Hessian.
#'
#'The mediator and outcome models (and exposure model for \code{type = "zm"} or \code{"zy"}) should be fitted using \code{glm} and can be of two types, probit models (\code{family = binomial(link = 'probit')})
#'for binary mediators or outcomes (exposures) and linear regression models (\code{family = gaussian}) for
#'continuous mediators or outcomes (exposures). Note that the exposure can either be binary or continuous, categorical exposures with more than two levels are not currently supported. The outcome model may contain exposure-mediator, exposure-covariate,
#'mediator-covariate and exposure-mediator-covariate interactions. The mediator model may contain exposure-covariate interactions.
#'All models may also contain interactions between covariates. Note, however that interactions may not be included in a model without
#'also including the main effects of the interacting variables. That is, interactions should be specified either as \code{X1*X2} or
#'\code{X1 + X2 + X1:X2}, not as \code{X1:X2} alone.
#'
#'To obtain results conditional on specific covariate values, these values should be provided through the \code{covariates} argument as a named list (see Examples).
#'The effects will be averaged over covariates not specified in the list.
#'
#'The total effect can be decomposed into a direct and indirect effect in different ways. Let z be the exposure value and z* the control (unexposed) value.
#'The default is to give the decomposition into the "pure direct effect" \code{E(Y(z,M(z*)))-E(Y(z*,M(z*)))} (here denoted NDE) and the "total indirect effect"
#'\code{E(Y(z,M(z)))-E(Y(z,M(z*)))} (denoted NIE). Setting \code{alt.decomposition=TRUE}
#'instead gives the decomposition into the "total direct effect" \code{E(Y(z,M(z)))-E(Y(z*,M(z)))} (here denoted NDE*) and "pure indirect effect"
#'\code{E(Y(z*,M(z)))-E(Y(z*,M(z*)))} (denoted NIE*).
#'
#'Standard errors for the effects are calculated using the delta method. Confidence intervals (CI) for (and p-values for tests of) the natural direct and indirect effects for each value of the
#'sensitivity parameter are constructed based on a normal approximation. Uncertainty intervals (UI) are constructed as the union of all CIs over the sensitivity parameter
#'vector.
#'
#'@return \code{sensmediation} returns an object of class \code{"effectsMed"}.
#'
#'The function \code{summary} (\code{\link{summary.effectsMed}}) gives a summary of the results in the form of a table with the estimated
#'effects and results of the sensitivity analysis. The function \code{plot} (\code{\link{plot.effectsMed}}) plots the estimated natural
#'indirect or direct effects with confidence intervals over the range of the sensitivity parameter.
#'
#'\item{call}{The matched call}
#'\item{Rho}{The sensitivity parameter vector.}
#'\item{type}{character, the type of confounding the sensitivity analysis is performed for.}
#'\item{coefs.sensmed}{a list with the output from \code{\link{coefs.sensmed}}}
#'\item{NIE}{matrix with the estimated NIEs (or NIE*s if \code{alt.decomposition=TRUE}) over the range of the sensitivity parameter \code{Rho}.}
#'\item{NDE}{matrix with the estimated NDEs (or NDE*s if \code{alt.decomposition=TRUE}) over the range of the sensitivity parameter \code{Rho}.}
#'\item{std.errs}{list with the standard errors of the NIE (NIE*), NDE (NDE*) and total effect over the range of the sensitivity parameter \code{Rho}.}
#'\item{CI}{a list with the confidence intervals of the NIE (NIE*), NDE (NDE*) and total effect over the range of the sensitivity parameter \code{Rho}.}
#'\item{UI}{matrix with the uncertainty intervals for the NIE (NIE*) and NDE (NDE*) over the range of the sensitivity parameter \code{Rho}.}
#'\item{conf.level}{numeric, the confidence level used for confidence intervals and uncertainty intervals.}
#'\item{covariates}{list of the covariate values that the effects are conditioned on.}
#'\item{exp.name}{character vector containing the name of the exposure variable.}
#'\item{med.name}{character vector containing the name of the mediator variable.}
#'\item{exp.value}{value of the exposure variable used as the exposure condition.}
#'\item{control.value}{value of the exposure variable used as the control (unexposed) condition.}
#'\item{alt.decomposition}{logical, indicating whether the alternative definitions of the direct and indirect effects have been used}
#'\item{med.model}{the mediator model input.}
#'\item{out.model}{the outcome model input.}
#'\item{betas}{list of the estimated mediator model parameters over \code{Rho}, with
#'\itemize{
#'\item \code{beta0} Intercept
#'\item \code{beta1} Exposure
#'\item \code{beta2} Covariates
#'\item \code{beta3} Exposure-covariate interactions
#'}
#'Components that are not included in the input mediator model are set to 0.}
#'\item{thetas}{list of the estimated outcome model parameters over \code{Rho}, with
#'\itemize{
#'\item \code{theta0} Intercept
#'\item \code{theta1} Exposure
#'\item \code{theta2} Mediator
#'\item \code{theta3} Exposure-mediator interaction
#'\item \code{theta4} Covariates
#'\item \code{theta5} Exposure-covariate interactions
#'\item \code{theta6} Mediator-covariate interactions
#'\item \code{theta7} Exposure-mediator-covariate interactions
#'}
#'Components that are not included in the input outcome model are set to 0.}
#'\item{part.deriv}{List with the partial derivatives of the NDE (Lambda), NIE (Gamma) and TE (Eta) wrt the mediator and outcome model parameters for each value of \code{Rho}. See \code{\link{partdevs}}.}
#'\item{sigma.thetabeta}{a list with the joint covariance matrix of the outcome and mediator model parameters for each value of \code{Rho}. Note that the covariance matrix is constructed for all estimated parameters listed in \code{betas} and \code{thetas} but that components not included in the input mediator and outcome models are set to 0.}
#'@export
#'@author Anita Lindmark
#'@references Lindmark, A., de Luna, X., Eriksson, M. (2018) Sensitivity Analysis for Unobserved Confounding of Direct and Indirect Effects Using Uncertainty Intervals, \emph{Statistics in Medicine}, \bold{37(10)}, pp 1744--1762, doi:10.1002/sim.7620.
#'
#' Lindmark A (2022). Sensitivity analysis for unobserved confounding in causal mediation analysis allowing for effect modification, censoring and truncation. \emph{Statistical Methods & Applications}, \bold{31}, pp 785--814, doi:10.1007/s10260-021-00611-4.
#'@seealso \code{\link{more.effects}} which can be used to calculate additional direct and indirect effects with sensitivity analysis using the same sensitivity parameter without running the optimization again.
#'@examples
#'
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
#' # Estimation of NIE, NDE and sensitivity analyses to mediator-outcome confounding:
#' # (note that the name of the exposure is "AF1" to match the name in coef(out.model))
#' sensmed <- sensmediation(m.model, o.model, exp.name = "AF1", med.name = "lowered.consc",
#'    Rho = c(0, 0.1))
#' summary(sensmed)
#' plot(sensmed)
#' plot(sensmed, effect = "direct")
#'
#' \dontrun{
#' # Conditional effects and sensitivity analysis to mediator-outcome confounding using
#' # more.effects():
#' sensmed.cond <- more.effects(sensmed.object = sensmed,
#'    covariates = list(sex = 1, age.cat = "70-79"))
#' summary(sensmed.cond)
#' }
#'
#' \dontrun{
#' ## Sensitivity analysis to exposure-mediator confounding:
#'   e.model <- glm(AF ~ age.cat + sex, data = RSdata,
#'      family = binomial(link = 'probit'))
#'
#'   sensmed.zm <- sensmediation(med.model = m.model, out.model = o.model,
#'      exp.model = e.model, type = "zm", Rho = seq(0, 0.5, 0.1), exp.name = "AF1",
#'      med.name = "lowered.consc")
#'
#'   summary(sensmed.zm)
#' }
#'
#' \dontrun{
#' # Additional effects using more.effects:
#' # Results with conf.level = 0.99:
#' sensmed.zm.99 <- more.effects(sensmed.object = sensmed.zm, conf.level = 0.99)
#' summary(sensmed.zm.99)
#' }
#'
#' \dontrun{
#' # Examples with simulated data, continuous exposure:
#'
#' require(mvtnorm)
#'
#' n <- 1000
#' set.seed(102677)
#'
#' x <- rnorm(n)
#' z <- -0.5 + 0.1*x + rnorm(n)
#' R <- 0.5
#' Sigma <- cbind(c(1,R), c(R,1))
#' epsilon <- rmvnorm(n, sigma = Sigma)
#' m <- -1.2 + 0.8*z + 0.13*x + epsilon[,1]
#' y <- -1 + 0.05*z + 3*m + 0.5*x + epsilon[,2]
#'
#' # Models:
#' z.model <- glm(z ~ x)
#' m.model2 <- glm(m ~ z + x)
#' y.model <- glm(y ~ z + m + x)
#'
#' ## Estimation of NIE, NDE. Note that the exposure condition is 2
#' ## so effects are calculated for a 2 unit increase of the exposure:
#' eff.contz <- sensmediation(med.model = m.model2, out.model = y.model,
#'            exp.name = "z", med.name = "m", control.value = 0, exp.value = 2)
#' summary(eff.contz)
#'}
#'


sensmediation <- function(med.model, out.model, exp.model = NULL, exp.name = NULL, med.name = NULL,
                           type = "my", Rho = 0, progress = TRUE, conf.level = 0.95, covariates = NULL,
                           alt.decomposition = FALSE, control.value = 0, exp.value = 1, covariance = NULL,
                           med.full = NULL, out.full = NULL, all.interactions = NULL, ...){

  cll <- match.call()

  effects <- list()

  # Checks to make sure that the exposure and mediator names have been input -----------------------------------
  if(is.null(exp.name))
    stop("The name of the exposure variable needs to be specified through the argument exp.name.")
  if(is.null(med.name))
    stop("The name of the mediator variable needs to be specified through the argument med.name.")

  # Checks to make sure that the models input are of the correct type ------------------------------------------
  if(class(med.model)[1] != "glm" | class(out.model)[1] != "glm")
    stop("All in-models need to be of class glm.")

  if(med.model$family$link != "probit" & med.model$family$family != "gaussian")
    stop("The combination of mediator and outcome model classes is not implemented.")
  if(out.model$family$link != "probit" & out.model$family$family != "gaussian")
    stop("The combination of mediator and outcome model classes is not implemented.")

  if(length(med.model$y) != length(out.model$y))
    stop("The number of observations in the mediator and outcome models do not match.")

  if(type != "zm" & type != "my" & type != "zy")
    stop("Invalid type, type should be one of 'my', 'zm' or 'zy'.")

  if(!is.null(covariance))
    warning("The covariance argument has been deprecated and is ignored (see the help file for the sensmediation function)")
  if(!is.null(med.full))
    warning("The med.full argument has been deprecated and is ignored (see the help file for the sensmediation function)")
  if(!is.null(out.full))
    warning("The out.full argument has been deprecated and is ignored (see the help file for the sensmediation function)")
  if(!is.null(all.interactions))
    warning("The all.interactions argument has been deprecated and is ignored (see the help file for the sensmediation function)")


  # Results for type="zm" --------------------------------------------------------------------------------------
  if(type == "zm"){
    # Checks of the exposure model
    if(is.null(exp.model))
      stop("An exposure model needs to be specified, see argument exp.model.")
    if(!is.null(exp.model)){
      if(class(exp.model)[1] != "glm")
        stop("exp.model needs to be of class glm.")
      if(exp.model$family$link != "probit" & exp.model$family$family != "gaussian")
        stop("The exposure model family should be binomial(probit) or gaussian.")
    }
    if(length(med.model$y) != length(exp.model$y))
      stop("The number of observations in the exposure and mediator models do not match.")

    # Maximum likelihood estimation of the regression parameters in the mediator (and exposure) model under
    # exposure-mediator confounding
    ML.object  <- coefs.sensmed(model.expl = exp.model, model.resp = med.model, Rho, progress = progress, ...)

    # Estimation of natural direct and indirect effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "zm", out.model = out.model, covariates = covariates, alt.decomposition = alt.decomposition,
                             exp.value = exp.value, control.value = control.value, exp.name = exp.name, med.name = med.name)


  }

  # Results for type="my" ----------------------------------------------------------------------------------------
  if(type == "my"){
    # Maximum likelihood estimation of the regression parameters in the mediator and outcome models under
    # mediator-outcome confounding
    ML.object  <- coefs.sensmed(model.expl = med.model, model.resp = out.model, Rho, progress = progress, ...)

    # Estimation of natural direct and indirect effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "my", exp.name = exp.name, med.name = med.name, covariates = covariates,
                              alt.decomposition = alt.decomposition, exp.value = exp.value, control.value = control.value)

    }

  # Results for type="zy" ----------------------------------------------------------------------------------------
  if(type == "zy"){
    # Checks of the exposure model
    if(is.null(exp.model))
      stop("An exposure model needs to be specified, see argument exp.model.")
    if(!is.null(exp.model)){
      if(class(exp.model)[1] != "glm")
        stop("exp.model needs to be of class glm.")
      if(exp.model$family$link != "probit" & exp.model$family$family != "gaussian")
        stop("The exposure model family should be binomial(probit) or gaussian.")
    }
    if(length(med.model$y) != length(exp.model$y))
      stop("The number of observations in the exposure and outcome models do not match.")

    # Maximum likelihood estimation of the regression parameters in the outcome (and exposure) model under
    # exposure-outcome confounding
    ML.object  <- coefs.sensmed(model.expl = exp.model, model.resp = out.model, Rho, progress = progress, ...)

    # Estimation of natural direct and indirect effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "zy", med.model = med.model, covariates = covariates, alt.decomposition = alt.decomposition,
                           exp.value = exp.value, control.value = control.value, exp.name = exp.name, med.name = med.name)

  }

  # Storage of components ---------------------------------------------------------------------------
  effects$out.model <- out.model
  effects$med.model <- med.model
  effects$NIE <- eff.calc$effects$NIE
  effects$NDE <- eff.calc$effects$NDE
  effects$std.errs <- eff.calc$std.errs
  effects$betas <- eff.calc$betas
  effects$thetas <- eff.calc$thetas
  effects$part.deriv <- eff.calc$part.deriv
  effects$sigma.thetabeta <- eff.calc$sigma.thetabeta
  effects$coefs.sensmed <- ML.object

  # Calculating confidence intervals for the natural direct, indirect and total effects ----------------------------
  effects$CI$CI.nie <- effects$CI$CI.nde <- effects$CI$CI.te <- matrix(ncol = 2, nrow = length(ML.object$Rho))
  dimnames(effects$CI$CI.nie) <- dimnames(effects$CI$CI.nde) <- dimnames(effects$CI$CI.te) <- list(paste(ML.object$Rho), c("lower", "upper"))

  effects$CI$CI.nie[,1] <- effects$NIE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nie
  effects$CI$CI.nie[,2] <- effects$NIE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nie
  effects$CI$CI.nde[,1] <- effects$NDE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nde
  effects$CI$CI.nde[,2] <- effects$NDE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nde
  effects$CI$CI.te[,1] <- effects$NIE + effects$NDE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.te
  effects$CI$CI.te[,2] <- effects$NIE + effects$NDE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.te

  # Calculating uncertainty intervals for the natural direct and indirect effect -----------------------------------
  effects$UI <- matrix(nrow = 2, ncol = 2)
  effects$UI[1,1] <- min(effects$CI$CI.nie[,1])
  effects$UI[1,2] <- max(effects$CI$CI.nie[,2])
  effects$UI[2,1] <- min(effects$CI$CI.nde[,1])
  effects$UI[2,2] <- max(effects$CI$CI.nde[,2])
  colnames(effects$UI) <- c("Lower", "Upper")
  if(alt.decomposition == TRUE)
    rownames(effects$UI) <- c("NIE*", "NDE*")
  else
    rownames(effects$UI) <- c("NIE", "NDE")

  # Storing relevant additional components -------------------------------------------------------------------------
  effects$covariates <- eff.calc$covariates
  effects$conf.level <- conf.level
  effects$Rho <- ML.object$Rho
  effects$type <- type
  effects$alt.decomposition <- alt.decomposition
  effects$med.name <- med.name
  effects$exp.name <- exp.name
  effects$control.value <- control.value
  effects$exp.value <- exp.value

  class(effects) <- "effectsMed"

  effects$call <- cll

  return(effects)

}


#'Estimate additional natural direct and indirect effects based on an \code{"effectsMed"} object
#'
#'Takes an \code{"effectsMed"} object and estimates additional natural direct and indirect effects, with a sensitivity analysis using
#'the same sensitivity parameter as in the original analysis, without having to redo the optimization to find the estimated regression coefficients.
#'The effects to be estimated are regulated through the arguments \code{covariates}, \code{alt.decomposition}, \code{exp.value} and \code{control.value}
#'as described in the documentation for \code{\link{sensmediation}}. The confidence level used is regulated through the argument \code{conf.level}.
#'
#'@param sensmed.object an object of class "effectsMed" for which additional effects are to be calculated.
#'@param conf.level the confidence level to be used for confidence intervals and uncertainty intervals.
#'@param covariates if conditional effects are to be estimated the list of covariate values (see \code{\link{sensmediation}}). Covariates not specified are marginalized over.
#'@param alt.decomposition logical indicating whether alternative definitions of the direct and indirect effects should be used (see \code{\link{sensmediation}}).
#'@param exp.value value of the exposure variable used as the exposure condition, default is to take the value stored in \code{sensmed.object}.
#'@param control.value value of the exposure variable used as the control (unexposed) condition, default is to take the value stored in \code{sensmed.object}.
#'@return \code{more.effects} returns an object of class \code{"effectsMed"}, see the documentation for \code{\link{sensmediation}}
#'for information.
#'
#'@export
#'@author Anita Lindmark
#'@seealso \code{\link{sensmediation}}
#'@examples
#'\dontrun{
#'
#' # Example with data from Riksstroke (the Swedish stroke register)
#'
#' data(RSdata)
#'
#' # Probit mediator and outcome models:
#' med.model <- glm(lowered.consc ~ AF + age.cat + sex, data = RSdata,
#' family = binomial(link = 'probit'))
#' out.model <- glm(cf.3mo ~ AF + lowered.consc + age.cat + sex, data = RSdata,
#' family = binomial(link = 'probit'))
#'
#' # First we estimate marginal NIE, NDE with sensitivity analyses to mediator-outcome
#' # confounding:
#' sensmed <- sensmediation(med.model, out.model, exp.name = "AF1", med.name = "lowered.consc",
#' Rho = seq(0, 0.5, 0.1))
#'
#' # Then we also estimate NIE, NDE conditional on male sex without reestimating the regression
#' # coefficients:
#' sensmed.cond <- more.effects(sensmed.object = sensmed, covariates = list(sex = 1))
#' summary(sensmed.cond)
#' plot(sensmed.cond)
#'}
#'


more.effects <- function(sensmed.object, conf.level = 0.95, covariates = NULL, alt.decomposition = FALSE, exp.value = NULL, control.value = NULL){

  cll <- match.call()

  effects <- sensmed.object

  type <- sensmed.object$type
  ML.object <- sensmed.object$coefs.sensmed

  if(is.null(control.value))
    control.value <- sensmed.object$control.value
  if(is.null(exp.value))
    exp.value <- sensmed.object$exp.value

  # Additional results for type="zm" ---------------------------------------------------------
  if(type == "zm"){

    # Estimating and storing the effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "zm", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, out.model = sensmed.object$out.model,
                             covariates = covariates, alt.decomposition = alt.decomposition,
                             exp.value = exp.value, control.value = control.value)



  }

  # Additional results for type="my" ---------------------------------------------------------
  if(type == "my"){
    # Estimating and storing the effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "my", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, covariates = covariates,
                             alt.decomposition = alt.decomposition, exp.value = exp.value,
                             control.value = control.value)
  }



  # Additional results for type="zy" ---------------------------------------------------------
  if(type == "zy"){
    # Estimating and storing the effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "zy", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, med.model = sensmed.object$med.model,
                             covariates = covariates, alt.decomposition = alt.decomposition,
                             exp.value = exp.value, control.value = control.value)

  }

  effects$NIE <- eff.calc$effects$NIE
  effects$NDE <- eff.calc$effects$NDE
  effects$std.errs <- eff.calc$std.errs
  effects$part.deriv <- eff.calc$part.deriv

  # Calculating confidence intervals for the natural direct and indirect effect and the total effect:
  effects$CI$CI.nie <- effects$CI$CI.nde <- effects$CI$CI.te <- matrix(ncol=2, nrow=length(ML.object$Rho))
  dimnames(effects$CI$CI.nie) <- dimnames(effects$CI$CI.nde ) <- dimnames(effects$CI$CI.te) <- list(paste(ML.object$Rho), c("lower", "upper"))

  effects$CI$CI.nie[,1] <- effects$NIE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nie
  effects$CI$CI.nie[,2] <- effects$NIE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nie
  effects$CI$CI.nde[,1] <- effects$NDE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nde
  effects$CI$CI.nde[,2] <- effects$NDE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.nde
  effects$CI$CI.te[,1] <- effects$NIE + effects$NDE - stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.te
  effects$CI$CI.te[,2] <- effects$NIE + effects$NDE + stats::qnorm(1 - (1 - conf.level)/2)*effects$std.errs$se.te

  # Calculating confidence intervals for the natural direct and indirect effect
  effects$UI <- matrix(nrow = 2, ncol = 2)
  effects$UI[1,1] <- min(effects$CI$CI.nie[,1])
  effects$UI[1,2] <- max(effects$CI$CI.nie[,2])
  effects$UI[2,1] <- min(effects$CI$CI.nde[,1])
  effects$UI[2,2] <- max(effects$CI$CI.nde[,2])
  colnames(effects$UI) <- c("Lower", "Upper")
  if(alt.decomposition==TRUE)
    rownames(effects$UI) <- c("NIE*", "NDE*")
  else
    rownames(effects$UI) <- c("NIE", "NDE")

  # Storing relevant additional components
  effects$covariates <- eff.calc$covariates
  effects$conf.level <- conf.level
  effects$alt.decomposition <- alt.decomposition

  effects$control.value <- control.value
  effects$exp.value <- exp.value


  class(effects) <- "effectsMed"

  effects$call <- cll

  return(effects)


}

