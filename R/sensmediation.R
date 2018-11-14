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
#'@param exp.name A character string indicating the name of the exposure variable used in the models.
#'@param med.name A character string indicating the name of the mediator used in the models.
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
#'\code{\link{maxLik}}, the default is to use the Newton-Raphson method and an analytic gradient and Hessian.
#'
#'The mediator and outcome models (and exposure model for \code{type = "zm"} or \code{"zy"}) should be fitted using \code{glm} and can be of two types, probit models (\code{family = binomial(link = 'probit')})
#'for binary mediators or outcomes (exposures) and linear regression models (\code{family = gaussian}) for
#'continuous mediators or outcomes (exposures). The outcome model may contain exposure-mediator, exposure-covariate,
#'mediator-covariate and exposure-mediator-covariate interactions. The mediator model may contain exposure-covariate interactions.
#'All models may also contain interactions between covariates. Note, however that interactions may not be included in a model without
#'also including the main effects of the interacting variables. That is, interactions should be specified either as \code{X1*X2} or
#'\code{X1 + X2 + X1:X2}, not as \code{X1:X2} alone.
#'
#'To obtain results conditional on specific covariate values, these values should be provided through the \code{covariates} argument as a named list. The names of the list
#'item(s) should match names from \code{names(med.model$coefficients)} and/or \code{names(out.model$coefficients)} (see Examples). The effects will be averaged
#'over covariates not specified in the list.
#'
#'The total effect can be decomposed into a direct and indirect effect in different ways. Let z be the exposure value and z* the control (unexposed) value.
#'The default is to give the decomposition into the "pure direct effect" \code{E(Y(z,M(z*)))-E(Y(z*,M(z*)))} (here denoted NDE) and the "total indirect effect"
#'\code{E(Y(z,M(z)))-E(Y(z,M(z*)))} (denoted NIE). Setting \code{alt.decomposition=TRUE}
#'instead gives the decomposition into the "total direct effect" \code{E(Y(z,M(z)))-E(Y(z*,M(z)))} (here denoted NDE*) and "pure indirect effect"
#'\code{E(Y(z*,M(z)))-E(Y(z*,M(z*)))} (denoted NIE*).
#'
#'Standard errors for the effects are calculated using the delta method.
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
#'\item{alt.decomposition}{logical, indicating whether the alternative definitions of the direct and indirect effects have been used}
#'\item{med.model}{if \code{type="zy"}, the mediator model input.}
#'\item{out.model}{if \code{type="zm"}, the outcome model input.}
#'@export
#'@author Anita Lindmark
#'@references Lindmark, A., de Luna, X., Eriksson, M. (2018) Sensitivity Analysis for Unobserved Confounding of Direct and Indirect Effects Using Uncertainty Intervals, \emph{Statistics in Medicine}, \bold{37(10)}, pp 1744--1762.
#'@seealso \code{\link{more.effects}} which can be used to calculate additional direct and indirect effects with sensitivity analysis using the same sensitivity parameter without running the optimization again.
#'@examples
#'
#' #Examples with simulated data
#'
#' #######################################
#' # Binary mediator, continuous outcome #
#' #######################################
#'
#' ###Binary exposure:
#'
#' ##Simulated data:
#' require(mvtnorm)
#' require(maxLik)
#'
#' n <- 1000
#' set.seed(102677)
#'
#' x <- rnorm(n)
#'
#' z.star <- -0.5 + 0.1*x + rnorm(n)
#' z <- ifelse(z.star > 0, 1, 0)
#'
#' R <- 0.5
#' Sigma <- cbind(c(1,R), c(R,1))
#' epsilon <- rmvnorm(n, sigma = Sigma)
#'
#' m.star <- -1.2 + 0.8*z + 0.13*x + epsilon[,1]
#' m <- ifelse(m.star > 0, 1, 0)
#'
#' y <- -1 + 0.05*z + 3*m + 0.5*x + epsilon[,2]
#'
#' #Models:
#' z.model <- glm(z ~ x, family = binomial(link = 'probit'))
#' m.model <- glm(m ~ z + x, family = binomial(link = 'probit'))
#' y.model <- glm(y ~ z + m + x)
#'
#'
#' ##Estimation of NIE, NDE and sensitivity analyses to mediator-outcome confounding:
#' effects.my <- sensmediation(med.model = m.model, out.model = y.model, exp.name = "z",
#'               med.name = "m", Rho = seq(0, 0.5, 0.1))
#' summary(effects.my)
#' summary(effects.my, non.sign = TRUE)
#' plot(effects.my)
#' plot(effects.my, effect="direct")
#'
#' ##Estimation of NIE, NDE and sensitivity analyses to exposure-mediator confounding:
#' \dontrun{
#'   effects.zm <- sensmediation(med.model = m.model, out.model = y.model, exp.model = z.model,
#'             type = "zm", Rho = seq(0, 0.5, 0.1), exp.name = "z", med.name = "m")
#'   summary(effects.zm)
#' }
#'
#' ##Additional effects using more.effects:
#' #Results with conf.level = 0.99:
#' effects.my.99 <- more.effects(sensmed.object = effects.my, conf.level = 0.99)
#' summary(effects.my.99)
#' #Conditional effects and sensitivity analysis to mediator-outcome confounding:
#' eff.my.cond <- more.effects(sensmed.object = effects.my, covariates = list(x = 1))
#' summary(eff.my.cond)
#'
#'
#' ###Continuous exposure:
#' #Models:
#' z.model.cont <- glm(z.star ~ x)
#' m.model.contz <- glm(m ~ z.star + x, family=binomial(link='probit'))
#' y.model.contz <- glm(y ~ z.star + m + x)
#'
#' ##Estimation of NIE, NDE and sensitivity analyses to mediator-outcome confounding:
#' eff.my.contz <- sensmediation(med.model = m.model.contz, out.model = y.model.contz,
#'            Rho = seq(0, 0.5, 0.1), exp.name = "z.star", med.name = "m",
#'            control.value = 0, exp.value = 2)
#' summary(eff.my.contz)
#'
#' ##Estimation of NIE, NDE and sensitivity analyses to exposure-mediator confounding:
#' eff.zm.contz <- sensmediation(med.model = m.model.contz, out.model = y.model.contz,
#'                  exp.model = z.model.cont, type = "zm", Rho = seq(0, 0.5, 0.1),
#'                  exp.name = "z.star", med.name = "m", control.value = 0, exp.value = 2)
#' summary(eff.zm.contz)
#'
#'

sensmediation <- function(med.model, out.model, exp.model = NULL, exp.name = NULL, med.name = NULL,
                           type = "my", Rho = 0, progress = TRUE, conf.level = 0.95, covariates = NULL,
                           alt.decomposition = FALSE, control.value = 0, exp.value = 1, covariance = NULL,
                           med.full = NULL, out.full = NULL, all.interactions = NULL, ...){

  cll <- match.call()

  effects <- list()

  # Checks to make sure that the necessary models have been input ----------------------------------------------
  if(is.null(exp.name))
    stop("The name of the exposure variable needs to be specified through the argument exp.model.")
  if(is.null(med.name))
    stop("The name of the mediator variable needs to be specified through the argument med.model.")

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

    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs
    effects$out.model <- out.model

  }

  # Results for type="my" ----------------------------------------------------------------------------------------
  if(type == "my"){
    # Maximum likelihood estimation of the regression parameters in the mediator and outcome models under
    # mediator-outcome confounding
    ML.object  <- coefs.sensmed(model.expl = med.model, model.resp = out.model, Rho, progress = progress, ...)

    # Estimation of natural direct and indirect effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "my", exp.name = exp.name, med.name = med.name, covariates = covariates,
                              alt.decomposition = alt.decomposition, exp.value = exp.value, control.value = control.value)

    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs

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

    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs
    effects$med.model <- med.model

  }

  # Storing the results of coefs.sensmed ---------------------------------------------------------------------------
  effects$coefs.sensmed <- ML.object

  # Calculating confidence intervals for the natural direct, indirect and total effects ----------------------------
  effects$CI$CI.nie <- matrix(ncol = 2, nrow = length(ML.object$Rho))
  dimnames(effects$CI$CI.nie) <- list(paste(ML.object$Rho), c("lower", "upper"))
  effects$CI$CI.nde <- effects$CI$CI.nie
  effects$CI$CI.te <- effects$CI$CI.nie

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
  effects$covariates <- covariates
  effects$conf.level <- conf.level
  effects$Rho <- ML.object$Rho
  effects$type <- type
  effects$alt.decomposition <- alt.decomposition
  effects$med.name <- med.name
  effects$exp.name <- exp.name

  class(effects) <- "effectsMed"

  effects$call <- cll

  return(effects)

}


#'Estimate additional natural direct and indirect effects based on an object from \code{sensmediation}
#'
#'Takes an object from \code{\link{sensmediation}} and estimates additional natural direct and indirect effects, with a sensitivity analysis using
#'the same sensitivity parameter as in the original analysis, without having to redo the optimization to find the estimated regression coefficients.
#'The effects to be estimated are regulated through the arguments \code{covariates} and \code{alt.decomposition} as described in the documentation
#'for \code{\link{sensmediation}}. The confidence level used is regulated through the argument \code{conf.level}.
#'
#'@param sensmed.object Object from \code{\link{sensmediation}} for which additional effects are to be calculated.
#'@param conf.level the confidence level to be used for confidence intervals and uncertainty intervals.
#'@param covariates if conditional effects are to be estimated the list of covariate values (see \code{\link{sensmediation}}). Covariates not specified are marginalized over.
#'@param alt.decomposition logical indicating whether alternative definitions of the direct and indirect effects should be used (see \code{\link{sensmediation}}).
#'@param exp.value value of the exposure variable used as the exposure condition, default is 1.
#'@param control.value value of the exposure variable used as the control (unexposed) condition, default is 0.
#'@return \code{more.effects} returns an object of class \code{"effectsMed"}.
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
#'\item{alt.decomposition}{logical, indicating whether the alternative definitions of the direct and indirect effects have been used}
#'@export
#'@author Anita Lindmark
#'@seealso \code{\link{sensmediation}}
#'@examples
#'\dontrun{
#'#First we estimate marginal NIE, NDE with sensitivity analyses to mediator-outcome confounding:
#'effects.my <- sensmediation(med.model = m.model, out.model = y.model, exp.name = "z",
#'                      med.name = "m", Rho = seq(0, 0.5, 0.1))
#'
#'#Then we want to do the same for conditional NIE, NDE without reestimating the regression
#'#coefficients:
#'effects.my.cond <- more.effects(sensmed.object = effects.my, covariates = list(x = 1))
#'summary(effects.my.cond)
#'plot(effects.my.cond)
#'}
#'


more.effects <- function(sensmed.object, conf.level = 0.95, covariates = NULL, alt.decomposition = FALSE, exp.value = 1, control.value = 0){
  
  cll <- match.call()
  
  effects <- list()
  
  type <- sensmed.object$type
  ML.object <- sensmed.object$coefs.sensmed
  
  # Additional results for type="zm" ---------------------------------------------------------
  if(type == "zm"){
    
    # Estimating and storing the effects and their standard errors
    out.model <- sensmed.object$out.model
    eff.calc <- calc.effects(ML.object, type = "zm", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, out.model = out.model, covariates = covariates,
                             alt.decomposition = alt.decomposition, exp.value = exp.value,
                             control.value = control.value)
    
    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs
    effects$out.model <- out.model
    
    
  }
  
  # Additional results for type="my" ---------------------------------------------------------
  if(type == "my"){
    # Estimating and storing the effects and their standard errors
    eff.calc <- calc.effects(ML.object, type = "my", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, covariates = covariates,
                             alt.decomposition = alt.decomposition, exp.value = exp.value,
                             control.value = control.value)
    
    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs
  }
  
  
  
  # Additional results for type="zy" ---------------------------------------------------------
  if(type == "zy"){
    # Estimating and storing the effects and their standard errors
    med.model <- sensmed.object$med.model
    eff.calc <- calc.effects(ML.object, type = "zy", exp.name = sensmed.object$exp.name,
                             med.name = sensmed.object$med.name, med.model = med.model,
                             covariates = covariates, alt.decomposition = alt.decomposition,
                             exp.value = exp.value, control.value = control.value)
    effects$NIE <- eff.calc$effects$NIE
    effects$NDE <- eff.calc$effects$NDE
    effects$std.errs <- eff.calc$std.errs
    effects$med.model <- med.model
  }
  
  
  
  effects$exp.name <- sensmed.object$exp.name
  effects$med.name <- sensmed.object$med.name
  
  # Storing the results of coefs.sensmed
  effects$coefs.sensmed <- ML.object
  
  # Calculating confidence intervals for the natural direct and indirect effect and the total effect:
  effects$CI$CI.nie <- matrix(ncol=2, nrow=length(ML.object$Rho))
  dimnames(effects$CI$CI.nie  ) <- list(paste(ML.object$Rho), c("lower", "upper"))
  effects$CI$CI.nde <- effects$CI$CI.nie
  effects$CI$CI.te <- effects$CI$CI.nie
  
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
  effects$covariates <- covariates
  effects$conf.level <- conf.level
  effects$Rho <- sensmed.object$Rho
  effects$type <- type
  effects$alt.decomposition <- alt.decomposition
  
  class(effects) <- "effectsMed"
  
  effects$call <- cll
  
  return(effects)
  
  
}

