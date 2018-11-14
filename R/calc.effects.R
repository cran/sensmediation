#'Function for estimation of natural direct and indirect effects and sensitivity analysis for unobserved mediator-outcome confounding
#'
#'Function to estimate natural direct and indirect effect estimates and standard errors (using the delta method) based on parametric regression models and perform sensitivity analysis for unobserved confounding.
#'Intended to be called through \code{\link{sensmediation}} (or \code{\link{more.effects}}),  not on its own.
#'
#'@param ML.object object from \code{\link{coefs.sensmed}}
#'@param type the type of confounding for which the sensitivity analysis is to be performed. \code{type = "my"},  the default, corresponds to unobserved mediator-outcome
#'confounding,  \code{type = "zm"} to exposure-mediator confounding and \code{type = "zy"} to exposure-outcome confounding.
#'@param exp.name A character string indicating the name of the exposure variable used in the models.
#'@param med.name A character string indicating the name of the mediator used in the models.
#'@param covariates if conditional effects are to be estimated the list of covariate values. Covariates not specified are marginalized over. For more information, see \code{\link{sensmediation}}.
#'@param alt.decomposition logical indicating whether alternative definitions of the direct and indirect effects should be used (for more information, see \code{\link{sensmediation}}).
#'@param exp.value value of the exposure variable used as the exposure condition, default is 1.
#'@param control.value value of the exposure variable used as the control (unexposed) condition, default is 0.
#'@param med.model If \code{type = "zy"},  fitted \code{\link{glm}} model object representing the mediator model at the basis of the estimation.
#'@param out.model If \code{type = "zm"},  fitted \code{\link{glm}} model object representing the outcome model at the basis of the estimation.
#'@return A list with elements:
#'\item{effects}{A list with elements \code{NIE} and \code{NDE},  row matrices with the estimated NIE and NDE (or NIE* and NDE* if \code{alt.decomposition = TRUE}) for each value of the sensitivity parameter \code{Rho}.}
#'\item{std.errs}{A list with elements \code{se.nie} and \code{se.nde},  row matrices with the estimated standard errors for the natural direct and indirect effects for the different values of the sensitivity parameter \code{Rho}.}
#'@author Anita Lindmark
#'@seealso \code{\link{sensmediation}}
#'@export


calc.effects <- function(ML.object, type="my", exp.name, med.name, covariates = NULL, alt.decomposition = FALSE,
                          exp.value = 1, control.value = 0, med.model = NULL, out.model = NULL)
{
  # Extracting the necessary components for estimation from the ML.object ---------------------------
  nrho <- length(ML.object$Rho)
  Rho <- ML.object$Rho

  if(type == "zm"){
    model.med <- ML.object$model.resp
    model.out <- out.model
  }
  if(type == "zy"){
    model.med <- med.model
    model.out <- ML.object$model.resp
  }
  if(type == "my"){
    model.med <- ML.object$model.expl
    model.out <- ML.object$model.resp
  }

  cb <- model.med$family$family == "gaussian" & model.out$family$link == "probit"
    # -----------------------------------------------------------------------------------------------

  # Finding the position of the different types of regression coefficients in the mediator and outcome models --
  pos.theta1 <- match(exp.name, names(model.out$coefficients))
  if(is.na(pos.theta1))
    stop("The exposure is either missing from the response model or incorrectly named in exp.name.")

  pos.beta1 <- match(exp.name, names(model.med$coefficients))
  if(is.na(pos.beta1))
    stop("The exposure is either missing from the mediator model or incorrectly named in exp.name.")

  pos.theta2 <- match(med.name, names(model.out$coefficients))
  if(is.na(pos.theta2))
    stop("The mediator is either missing from the response model or incorrectly named in med.name.")

  # Identifying positions with interaction terms.
  out.split <- strsplit(names(model.out$coefficients), ":") # Split model terms by ":"
  med.split <- strsplit(names(model.med$coefficients), ":")
  logical.int.out <- sapply(lapply(out.split, length), ">", 1)
  logical.int.med <- sapply(lapply(med.split, length), ">", 1)
  int.out.split <- out.split[logical.int.out]
  int.med.split <- med.split[logical.int.med]
  pos.out <- c(1:length(out.split))
  pos.med <- c(1:length(med.split))
  pos.int.out <- pos.out[logical.int.out]
  pos.int.med <- pos.med[logical.int.med]

  # Logical vectors, do the interactions involve the mediator, exposure or both:
  logical.mint.out <- unlist(lapply(int.out.split, FUN = "%in%" , x = med.name))
  logical.zint.out <- unlist(lapply(int.out.split, FUN = "%in%" , x = exp.name))
  logical.zmint.out <- logical.mint.out == TRUE & logical.zint.out == TRUE

  pos.beta3 <- pos.int.med[unlist(lapply(int.med.split, FUN = "%in%" , x = exp.name))]

  pos.zmint <- pos.int.out[logical.zmint.out]
  pos.theta3 <- pos.zmint[which(unlist(lapply(int.out.split[logical.zmint.out], length)) == 2)]
  pos.theta7 <- pos.zmint[which(unlist(lapply(int.out.split[logical.zmint.out], length)) > 2)]

  pos.theta6 <- pos.int.out[logical.mint.out]
  if(length(pos.zmint))
    pos.theta6 <- pos.theta6[-which(pos.theta6%in%pos.zmint)]

  pos.theta5 <- pos.int.out[logical.zint.out]
  if(length(pos.zmint))
    pos.theta5 <- pos.theta5[-which(pos.theta5%in%pos.zmint)]

  indic.interc.out <- ifelse(attributes(model.out$terms)$intercept == 0, 0, 1)
  pos.theta4 <- pos.out[-which(pos.out%in%c(indic.interc.out, pos.theta1, pos.theta2, pos.theta3, pos.theta5,
                                            pos.theta6, pos.theta7))]
  pos.theta4 <- sort(pos.theta4)

  indic.interc.med <- ifelse(attributes(model.med$terms)$intercept == 0, 0, 1)
  pos.beta2 <- pos.med[-which(pos.med%in%c(indic.interc.med, pos.beta1, pos.beta3))]
  pos.beta2 <- sort(pos.beta2)
  # -----------------------------------------------------------------------------------------------

  # If sensitivity analysis to unobserved exposure-mediator confounding ---------------------------
  if(type == "zm"){
    # Extracting mediator and outcome regression parameters
    d.medcoef <- nrow(ML.object$coef)
    d.outcoef <- length(out.model$coef)

    medcoefs <- ML.object$coef # The matrix of mediator model coefficients (not including any sigmas)
    outcoefs <- matrix(out.model$coef, nrow = d.outcoef, ncol = nrho) # The matrix of outcome model coefficients.
    colnames(outcoefs) <- paste(Rho)                                  # Repeats the glm-coefficients length(Rho) times.
    rownames(outcoefs) <- names(out.model$coef)

    # Extracting covariance matrices for the mediator and outcome regression parameters:
    sigma.out <- stats::vcov(out.model) # Covariance matrix for the outcome model

    if(cb){
      d.medcoef <- d.medcoef + 1
      sigma.eta <- ML.object$sigma.res.resp
    }

    ph <- matrix(0, nrow = d.medcoef + d.outcoef, ncol = d.medcoef + d.outcoef)
    ph[1:d.outcoef, 1:d.outcoef] <- sigma.out
    sigma.pars <- lapply(1:nrho, function(x) x <- ph)
    sigma.pars <- lapply(1:nrho, function(x){
      sigma.pars[[x]][(d.outcoef + 1):(d.outcoef + d.medcoef), (d.outcoef + 1):(d.outcoef + d.medcoef)] <- ML.object$sigmas[[x]][1:(d.medcoef), 1:(d.medcoef)]; sigma.pars[[x]]})

  }
  # -----------------------------------------------------------------------------------------------

  # If sensitivity analysis to unobserved exposure-outcome confounding ----------------------------
  if(type=="zy"){
    # Extracting mediator and outcome regression parameters
    d.medcoef <- length(med.model$coef)
    d.outcoef <- nrow(ML.object$coef)

    medcoefs <- matrix(med.model$coef, nrow = d.medcoef, ncol = nrho)
    colnames(medcoefs) <- paste(Rho)
    rownames(medcoefs) <- names(med.model$coef)
    outcoefs <- ML.object$coef

    # Extracting covariance matrices for the mediator and outcome regression parameters:
    sigma.med <- stats::vcov(med.model)

    if(model.med$family$family=="gaussian" & model.out$family$link=="probit"){
      ph <- matrix(0, nrow = d.medcoef + d.outcoef + 1, ncol = d.medcoef + d.outcoef + 1)
      ph[(d.outcoef + 1):(d.outcoef + d.medcoef), (d.outcoef + 1):(d.outcoef + d.medcoef)] <- sigma.med
      ph[(d.outcoef + d.medcoef + 1), (d.outcoef + d.medcoef + 1)] <- summary(med.model)$dispersion/(2*med.model$df.residual)
      sigma.eta <- matrix(sqrt(summary(med.model)$dispersion), nrow = nrho, ncol = 1)
    }
    else{
      ph <- matrix(0, nrow = d.medcoef + d.outcoef, ncol = d.medcoef + d.outcoef)
      ph[(d.outcoef + 1):(d.outcoef + d.medcoef), (d.outcoef + 1):(d.outcoef + d.medcoef)] <- sigma.med

      }

    sigma.pars <- lapply(1:nrho, function(x) x <- ph)
    sigma.pars <- lapply(1:nrho, function(x){
      sigma.pars[[x]][1:d.outcoef, 1:d.outcoef] <- ML.object$sigmas[[x]][1:d.outcoef, 1:d.outcoef]; sigma.pars[[x]]})

    }

  # -----------------------------------------------------------------------------------------------

  # If sensitivity analysis to unobserved mediator-outcome confounding ----------------------------
  if(type == "my"){
    # Extracting mediator and outcome regression parameters
    d.medcoef <- nrow(ML.object$expl.coef)
    d.outcoef <- nrow(ML.object$coef)

    medcoefs <- ML.object$expl.coef
    outcoefs <- ML.object$coef

    # Extracting covariance matrices for the mediator and outcome regression parameters:
    if(model.out$family$link == "probit")
      sigma.pars <- ML.object$sigmas

    if(cb)
      sigma.eta <- ML.object$sigma.res.expl

    if(model.med$family$link == "probit" & model.out$family$family == "gaussian")
      sigma.pars <- lapply(1:nrho, function(x) ML.object$sigmas[[x]][-(d.outcoef + 1), -(d.outcoef + 1)])

    if(model.med$family$family == "gaussian" & model.out$family$family == "gaussian")
      sigma.pars <- lapply(1:nrho, function(x) ML.object$sigmas[[x]][-c((d.outcoef + 1), (d.outcoef + d.medcoef + 2)), -c((d.outcoef + 1), (d.outcoef + d.medcoef + 2))])

  }
  # -----------------------------------------------------------------------------------------------

  ##### Intercepts #####
  beta0 <- ifelse(rep(attributes(model.med$terms)$intercept, nrho) == 0, rep(0, nrho),
                  medcoefs[1, ])
  names(beta0) <- paste(Rho)
  theta0 <- ifelse(rep(attributes(model.out$terms)$intercept, nrho) == 0, rep(0, nrho),
                  outcoefs[1, ])
  names(theta0) <- paste(Rho)

  ##### Exposure and mediator (main effects) #####
  theta1 <- outcoefs[pos.theta1, ]
  names(theta1) <- paste(Rho)
  theta2 <- outcoefs[pos.theta2, ]
  names(theta2) <- paste(Rho)
  beta1 <- medcoefs[pos.beta1, ]
  names(beta1) <- paste(Rho)

  ##### Covariate coefficients #####
  covariate.coefs <- function(model, pos, coefs, nrho){
    if(length(pos) == 0){
      cov.coefs <- matrix(0, nrow = 1, ncol = nrho)
      covars <- rep(0, length(model$y))
    }
    if(length(pos) > 1){
      if(nrho == 1){
        cov.coefs <- matrix(coefs[pos,])
        rownames(cov.coefs) <- names(model$coefficients[pos])
      }
      else{
        cov.coefs <- coefs[pos, ]
      }

      covars <- stats::model.matrix(model)[, pos]

    }

    if(length(pos) == 1){
      cov.coefs <- t(as.matrix(coefs[pos, ] ))
      rownames(cov.coefs) <- names(model$coefficients[pos])
      covars <- stats::model.matrix(model)[, pos]
    }
    return(list(cov.coefs, covars))
  }

  theta4 <- covariate.coefs(model.out, pos.theta4, outcoefs, nrho)[[1]]
  colnames(theta4) <- paste(Rho)
  beta2 <- covariate.coefs(model.med, pos.beta2, medcoefs, nrho)[[1]]
  colnames(beta2) <- paste(Rho)
  x.med <- as.matrix(covariate.coefs(model.med, pos.beta2, medcoefs, nrho)[[2]])
  x.out <- as.matrix(covariate.coefs(model.out, pos.theta4, outcoefs, nrho)[[2]])
  colnames(x.med) <- rownames(beta2)
  colnames(x.out) <- rownames(theta4)

  ##### ZM interaction #####
  theta3 <- ifelse(rep(length(pos.theta3), nrho) == 0, rep(0, nrho),
         outcoefs[pos.theta3, ])
  names(theta3) <- paste(Rho)

  ##### Interactions btw covariates and exposure and/or mediator #####
  interactions <- function(model, pos, cov.coef, exp.name, med.name, coefs){

    if(length(pos) > ncol(cov.coef))
      stop("Check the mediator and outcome models. Interactions are only allowed if the main effects are also included in the model.")

    int.coefs <- matrix(0, nrow = ncol(cov.coef), ncol = ncol(coefs))
    matches <- integer(0)

    if(length(pos) == ncol(cov.coef)){
      int.coefs <- coefs[pos, ]
      matches <- c(1:length(pos))
      if(ncol(cov.coef) == 1)
        int.coefs <- t(as.matrix((int.coefs)))
      if(ncol(cov.coef > 1 & ncol(coefs) == 1))
        int.coefs <- as.matrix((int.coefs))
    }

    if(length(pos) > 0 & length(pos) != ncol(cov.coef)){
      names.split <- unlist(strsplit(names(model$coefficients)[pos],
                                     paste(exp.name,":", sep="")))
      names.split <- unlist(strsplit(names.split, paste(":", exp.name, sep="")))
      names.split <- unlist(strsplit(names.split, paste(med.name,":", sep="")))
      names.split <- unlist(strsplit(names.split, paste(":", med.name, sep="")))
      matches <- which(!is.na(match(colnames(cov.coef), names.split)))
      int.coefs[matches, ] <- coefs[pos, ]
    }

    return(list("int.coefs" = int.coefs, "matches" = matches))
  }

  theta5 <- interactions(model.out, pos.theta5, x.out, exp.name, med.name, outcoefs)$int.coefs
  rownames(theta5) <- paste(exp.name, ":", rownames(theta4), sep = "")
  colnames(theta5) <- paste(Rho)
  matches.theta5 <- interactions(model.out, pos.theta5, x.out, exp.name, med.name, outcoefs)$matches

  theta6 <- interactions(model.out, pos.theta6, x.out, exp.name, med.name, outcoefs)$int.coefs
  rownames(theta6) <- paste(med.name, ":", rownames(theta4), sep = "")
  colnames(theta6) <- paste(Rho)
  matches.theta6 <- interactions(model.out, pos.theta6, x.out, exp.name, med.name, outcoefs)$matches

  theta7 <- interactions(model.out, pos.theta7, x.out, exp.name, med.name, outcoefs)$int.coefs
  rownames(theta7) <- paste(exp.name, ":", med.name, ":", rownames(theta4), sep = "")
  colnames(theta7) <- paste(Rho)
  matches.theta7 <- interactions(model.out, pos.theta7, x.out, exp.name, med.name, outcoefs)$matches

  beta3 <- interactions(model.med, pos.beta3, x.med, exp.name, med.name, medcoefs)$int.coefs
  rownames(beta3) <- paste(exp.name, ":", rownames(beta2), sep = "")
  colnames(beta3) <- paste(Rho)
  matches.beta3 <- interactions(model.med, pos.beta3, x.med, exp.name, med.name, medcoefs)$matches


  betas <- list("beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3)
  thetas <- list("theta0" = theta0, "theta1" = theta1, "theta2" = theta2, "theta3" = theta3,
                 "theta4" = theta4, "theta5" = theta5, "theta6" = theta6, "theta7" = theta7)
  # -----------------------------------------------------------------------------------------------

  # Reordering the covariance matrices and augmenting them with zeroes (if necessary) -----------------------------
  pos.reord <- c(indic.interc.out, pos.theta1, pos.theta2, pos.theta3, pos.theta4, pos.theta5,
                 pos.theta6, pos.theta7, indic.interc.med*(length(model.out$coefficients)+1),
                 c(pos.beta1, pos.beta2, pos.beta3) + length(model.out$coefficients),
                 ifelse(cb, length(model.out$coefficients) + length(model.med$coefficients) + 1, 0) )

  sigma.pars <- lapply(sigma.pars, function(x) x[pos.reord, pos.reord])

  dim.full.out <- 4 + 4*nrow(theta4)
  dim.full.med <- 2 + 2*nrow(beta2) + ifelse(cb, 1, 0)
  dim.full <- dim.full.out + dim.full.med

  if(length(pos.reord) != dim.full){
    ph.full <- matrix(0, nrow = dim.full, ncol = dim.full)
    sigma.full <- lapply(1:nrho, function(x) x <- ph.full)

    int.m <- indic.interc.med *( dim.full.out + 1)

    pos.sigma <- c(indic.interc.out, 2, 3, ifelse(length(pos.theta3)==0, 0, 4),
                   seq(1, length(pos.theta4), length.out = length(pos.theta4)) + 4,
                   length(pos.theta4) + 4 + matches.theta5, length(pos.theta4)*2 + 4 + matches.theta6,
                   length(pos.theta4)*3 + 4 + matches.theta7, int.m, dim.full.out + 2,
                   dim.full.out + 2 + seq(1, length(pos.beta2), length.out = length(pos.beta2)),
                   length(pos.beta2) + dim.full.out + 2 + matches.beta3, ifelse(cb, dim.full, 0)  )

    sigma.full <- lapply(1:nrho, function(x){ sigma.full[[x]][pos.sigma, pos.sigma] <- sigma.pars[[x]]; sigma.full[[x]] })

    sigma.pars <- sigma.full
  }
  # -----------------------------------------------------------------------------------------------

  # If conditional effects are to be calculated ---------------------------------------------------------------
  if(!is.null(covariates)){

    med.w <- match(colnames(x.med), names(covariates), nomatch = 0)
    w.med <- match(names(covariates)[med.w], rownames(beta2), nomatch = 0)

    out.w <- match(colnames(x.out), names(covariates), nomatch = 0)
    w.out <- match(names(covariates)[out.w], rownames(theta4), nomatch = 0)

    covariates <- unlist(covariates)

    if(sum(w.med) == 0 & sum(w.out) == 0)
      stop("The names listed in the covariates argument do not match any of the covariates in the mediator or outcome models.")

    if(sum(w.med) != 0)
      x.med[, w.med] <- matrix(rep(covariates[med.w], each = nrow(x.med)), ncol = length(covariates[med.w]))

    if(sum(w.out) != 0)
      x.out[, w.out] <- matrix(rep(covariates[out.w], each = nrow(x.out)), ncol = length(covariates[out.w]))

  }
  # -----------------------------------------------------------------------------------------------

  # Estimation of effects and standard errors depending on the types of mediator and outcome models -----------------
  if(model.med$family$link == "probit"){
    if(model.out$family$link == "probit"){
      effects <- eff.bb(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value, control.value)
      stderr <- stderr.bb(Rho, betas, thetas, sigma.pars, x.med, x.out, alt.decomposition, exp.value, control.value)

      }
    if(model.out$family$family == "gaussian"){
      effects <- eff.bc(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value, control.value)
      stderr <- stderr.bc(Rho, betas, thetas, sigma.pars, x.med, x.out,
                           alt.decomposition, exp.value, control.value)
    }
  }

  if(model.med$family$family == "gaussian"){
    if(model.out$family$link == "probit"){
      effects <- eff.cb(Rho, betas, thetas, sigma.eta, x.med, x.out, alt.decomposition, exp.value, control.value)
      stderr <- stderr.cb(Rho, betas, thetas, sigma.eta, sigma.pars, x.med, x.out, alt.decomposition,
                           exp.value, control.value)
    }
    if(model.out$family$family == "gaussian"){
      effects <- eff.cc(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value, control.value)
      stderr <- stderr.cc(Rho, betas, thetas, sigma.pars, x.med, x.out, alt.decomposition, exp.value, control.value)
    }
  }
  # -----------------------------------------------------------------------------------------------

  return(list("effects" = effects, "std.errs" = stderr))

}

#'Functions to calculate natural direct and indirect effects.
#'
#'Functions used to calculate natural direct and indirect effects based on the estimated regression parameters. Called by \code{\link{calc.effects}}.
#'The functions are named according to the convention \code{eff."mediator model type""outcome model type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression.
#'
#'@param Rho The sensitivity parameter vector.
#'@param betas List of mediator regression parameters
#'@param thetas List of outcome regression parameters
#'@param sigma.eta For a continuous mediator and binary outcome, matrix with the estimated residual standard deviation for the mediator model over the range of \code{Rho}.
#'@param x.med Mediator covariate matrix for which to calculate standard errors
#'@param x.out Outcome covariate matrix for which to calculate standard errors
#'@param alt.decomposition logical indicating whether or not alternative definitions of the direct and indirect effects should be used.
#'@param exp.value value of the exposure variable used as the exposure condition.
#'@param control.value value of the exposure variable used as the control (unexposed) condition.
#'@name effects
NULL

#'@rdname effects
#'@export
eff.bb <- function(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value,
                   control.value){
  nrho <- length(Rho)
  NIE <- matrix(nrow = 1, ncol = nrho)
  colnames(NIE) <- paste(Rho)
  NDE <- NIE

  rownames(NIE) <- ifelse(alt.decomposition == TRUE, c("NIE*"), c("NIE"))
  rownames(NDE) <- ifelse(alt.decomposition == TRUE, c("NDE*"), c("NDE"))

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:nrho){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th0 <- thetas$theta0[i]
    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th4 <- thetas$theta4[, i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    probs.med.ie <- stats::pnorm(b0 + b1*exp.value + x.med%*%(b2 + b3*exp.value))-
      stats::pnorm(b0 + b1*control.value + x.med%*%(b2 + b3*control.value))

    probs.med.de <- stats::pnorm(b0 + b1*t.de + x.med%*%(b2 + b3*t.de))

    probs.out.ie <-  stats::pnorm(th0 + th2 + (th1 + th3)*t.ie + x.out%*%(th4 + th5*t.ie +
                                                                     th6 + th7*t.ie))-stats::pnorm(th0 + th1*t.ie + x.out%*%(th4 + th5*t.ie))

    probs.out.de1 <- stats::pnorm(th0 + th1*exp.value + x.out%*%(th4 + exp.value*th5)) -
      stats::pnorm(th0 + th1*control.value + x.out%*%(th4 + control.value*th5))

    probs.out.de2 <- stats::pnorm(th0 + th2 + exp.value*(th1 + th3) + x.out%*%(th4 +
                                                                          exp.value*th5 + th6 + exp.value*th7)) -
      stats::pnorm(th0 + th2 + control.value*(th1 + th3) + x.out%*%(th4 + control.value*th5 +
                                                               th6 + control.value*th7))


    NIE[, i] <- mean(probs.out.ie*probs.med.ie)
    NDE[, i] <- mean(probs.out.de1*(1 - probs.med.de )  +  probs.out.de2*probs.med.de )


  }
  return(list("NIE" = NIE, "NDE" = NDE))
}

#'@rdname effects
#'@export
eff.bc <- function(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value, control.value){
  nrho <- length(Rho)
  NIE <- matrix(nrow = 1, ncol = nrho)
  colnames(NIE) <- paste(Rho)
  NDE <- NIE

  rownames(NIE) <- ifelse(alt.decomposition == TRUE, c("NIE*"), c("NIE"))
  rownames(NDE) <- ifelse(alt.decomposition == TRUE, c("NDE*"), c("NDE"))

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  diff <- exp.value - control.value
  for(i in 1:nrho){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    med.ie <- stats::pnorm(b0 + b1*exp.value + x.med%*%(b2  +  exp.value*b3))-
      stats::pnorm(b0 + b1*control.value + x.med%*%(b2  +  control.value*b3))
    med.de <- stats::pnorm(b0 + b1*t.de + x.med%*%(b2  +  t.de*b3))

    out.ie <- th2  +  th3*t.ie  +  x.out%*%(th6  +  t.ie*th7)
    out.de <- th3*diff  +  x.out%*%th7*diff

    NIE[, i] <- mean(out.ie*med.ie)
    NDE[, i] <- mean(th1*diff +  x.out%*%th5*diff  +  out.de*med.de )


  }
  return(list("NIE" = NIE, "NDE" = NDE))
}

#'@rdname effects
#'@export
eff.cb <- function(Rho, betas, thetas, sigma.eta, x.med, x.out, alt.decomposition,
                   exp.value, control.value){
  nrho <- length(Rho)
  NIE <- matrix(nrow = 1, ncol = nrho)
  colnames(NIE) <- paste(Rho)
  NDE <- NIE

  rownames(NIE) <- ifelse(alt.decomposition == TRUE, c("NIE*"), c("NIE"))
  rownames(NDE) <- ifelse(alt.decomposition == TRUE, c("NDE*"), c("NDE"))

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:nrho){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th0 <- thetas$theta0[i]
    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th4 <- thetas$theta4[, i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]
    sigma.e <- sigma.eta[i, ]

    med.ie.tr <- b0 + b1*exp.value + x.med%*%(b2  +  exp.value*b3)
    med.ie.cont <- b0 + b1*control.value + x.med%*%(b2  +  control.value*b3)
    med.de <- b0 + b1*t.de + x.med%*%(b2  +  t.de*b3)

    out.ie <- th2  +  th3*t.ie + x.out%*%(th6  +  t.ie*th7)
    out.de.tr <- th2  +  th3*exp.value + x.out%*%(th6  +  exp.value*th7)
    out.de.cont <- th2  +  th3*control.value + x.out%*%(th6  +  control.value*th7)

    denom.ie <- sqrt(sigma.e^2*out.ie^2 + 1)

    nie <- stats::pnorm((th0 + th1*t.ie + out.ie*med.ie.tr +  x.out%*%(th4 + th5*t.ie) )/denom.ie )  -
      stats::pnorm((th0 + th1*t.ie + out.ie*med.ie.cont +  x.out%*%(th4 + th5*t.ie) )/denom.ie )

    denom.de.tr <- sqrt(sigma.e^2*out.de.tr^2 + 1)
    denom.de.cont <- sqrt(sigma.e^2*out.de.cont^2 + 1)

    nde <- stats::pnorm((th0 + th1*exp.value + out.de.tr*med.de +  x.out%*%(th4 + th5*exp.value) )/denom.de.tr )  -
      stats::pnorm((th0 + th1*control.value + out.de.cont*med.de +  x.out%*%(th4 + th5*control.value) )/denom.de.cont )

    NIE[, i] <- mean(nie)
    NDE[, i] <- mean(nde)

  }
  return(list("NIE" = NIE, "NDE" = NDE))
}

#'@rdname effects
#'@export
eff.cc <- function(Rho, betas, thetas, x.med, x.out, alt.decomposition, exp.value, control.value){
  nrho <- length(Rho)
  NIE <- matrix(nrow = 1, ncol = nrho)
  colnames(NIE) <- paste(Rho)
  NDE <- NIE

  rownames(NIE) <- ifelse(alt.decomposition == TRUE, c("NIE*"), c("NIE"))
  rownames(NDE) <- ifelse(alt.decomposition == TRUE, c("NDE*"), c("NDE"))

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  diff <- exp.value - control.value
  for(i in 1:nrho){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    med.ie <- b1*diff + x.med%*%b3*diff
    med.de <- b0 + b1*t.de + x.med%*%(b2 + t.de*b3)

    out.ie <- th2 + th3*t.ie + x.out%*%(th6 + t.ie*th7)
    out.de <- th3*diff + x.out%*%th7*diff

    NIE[,i]<- mean(out.ie*med.ie)
    NDE[,i]<- mean(th1*diff+ x.out%*%th5*diff + out.de*med.de )

  }
  return(list("NIE"=NIE,"NDE"=NDE))
}

#'Functions to calculate standard errors of the direct, indirect and total effects using the delta method.
#'
#'Functions used to calculate standard errors of the direct, indirect and total effects using the delta method. Called by \code{\link{calc.effects}}.
#'The functions are named according to the convention \code{stderr."mediator model type""outcome model type"} where \code{b}
#'stands for binary probit regression and \code{c} stands for linear regression.
#'
#'@param Rho The sensitivity parameter vector.
#'@param betas List of mediator regression parameters
#'@param thetas List of outcome regression parameters
#'@param sigma.eta For a continuous mediator and binary outcome, matrix with the estimated residual standard deviation for the mediator model over the range of \code{Rho}.
#'@param sigma.pars List of covariance matrices for the mediator and outcome regression parameters
#'@param x.med Mediator covariate matrix for which to calculate standard errors
#'@param x.out Outcome covariate matrix for which to calculate standard errors
#'@param alt.decomposition logical indicating whether or not alternative definitions of the direct and indirect effects should be used.
#'@param exp.value value of the exposure variable used as the exposure condition.
#'@param control.value value of the exposure variable used as the control (unexposed) condition.
#'@name stderrs
NULL


#'@rdname stderrs
#'@export
stderr.bb <- function(Rho, betas, thetas, sigma.pars, x.med, x.out, alt.decomposition, exp.value, control.value){

  se.nie <- matrix(nrow = 1, ncol = length(Rho))
  colnames(se.nie) <- paste(Rho)
  se.nde <- se.nie
  se.te <- se.nie
  n <- nrow(x.out)

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:length(Rho)){

    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th0 <- thetas$theta0[i]
    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th4 <- thetas$theta4[, i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    part.derivs <- partdevs.bb(beta0 = b0, beta1 = b1, beta2 = b2, beta3 = b3, theta0 = th0, theta1 = th1,
                               theta2 = th2, theta3 = th3, theta4 = th4, theta5 = th5, theta6 = th6,
                               theta7 = th7, x.med, x.out, t.de, t.ie, exp.value, control.value)

    se.nie[, i] <- sqrt(t(part.derivs$Gamma)%*%sigma.pars[[i]]%*%part.derivs$Gamma)
    se.nde[, i] <- sqrt(t(part.derivs$Lambda)%*%sigma.pars[[i]]%*%part.derivs$Lambda)
    se.te[, i] <- sqrt(t(part.derivs$Eta)%*%sigma.pars[[i]]%*%part.derivs$Eta)
  }

  ses <- list("se.nie" = se.nie, "se.nde" = se.nde, "se.te" = se.te)
  return(ses)

}

#'@rdname stderrs
#'@export
stderr.bc <- function(Rho, betas, thetas, sigma.pars, x.med, x.out,
                      alt.decomposition, exp.value, control.value){

  se.nie <- matrix(nrow = 1, ncol = length(Rho))
  colnames(se.nie) <- paste(Rho)
  se.nde <- se.nie
  se.te <- se.nie

  diff <- exp.value - control.value

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:length(Rho)){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    part.derivs <- partdevs.bc(beta0 = b0, beta1 = b1, beta2 = b2, beta3 = b3, theta2 = th2, theta3 = th3,
                               theta6 = th6, theta7 = th7, x.med = x.med, x.out = x.out, t.de = t.de,
                               t.ie = t.ie, exp.value = exp.value, control.value = control.value)

    se.nie[, i] <- sqrt(t(part.derivs$Gamma)%*%sigma.pars[[i]]%*%part.derivs$Gamma)
    se.nde[, i] <- sqrt(t(part.derivs$Lambda)%*%sigma.pars[[i]]%*%part.derivs$Lambda)
    se.te[, i] <- sqrt(t(part.derivs$Eta)%*%sigma.pars[[i]]%*%part.derivs$Eta)
  }

  ses <- list("se.nie" = se.nie, "se.nde" = se.nde, "se.te" = se.te)
  return(ses)

}

#'@rdname stderrs
#'@export
stderr.cb <- function(Rho, betas, thetas, sigma.eta, sigma.pars, x.med, x.out, alt.decomposition,
                      exp.value, control.value){

  se.nie <- matrix(nrow = 1, ncol = length(Rho))
  colnames(se.nie) <- paste(Rho)
  se.nde <- se.nie
  se.te <- se.nie

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:length(Rho)){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th0 <- thetas$theta0[i]
    th1 <- thetas$theta1[i]
    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th4 <- thetas$theta4[, i]
    th5 <- thetas$theta5[, i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    part.derivs <- partdevs.cb(beta0 = b0, beta1 = b1, beta2 = b2, beta3 = b3, theta0 = th0, theta1 = th1,
                               theta2 = th2, theta3 = th3, theta4 = th4, theta5 = th5, theta6 = th6,
                               theta7 = th7, sigma.eta = sigma.eta[i, ], x.med = x.med, x.out = x.out,
                               t.de = t.de, t.ie = t.ie, exp.value = exp.value, control.value = control.value)

    se.nie[, i] <- sqrt(t(part.derivs$Gamma)%*%sigma.pars[[i]]%*%part.derivs$Gamma)
    se.nde[, i] <- sqrt(t(part.derivs$Lambda)%*%sigma.pars[[i]]%*%part.derivs$Lambda)
    se.te[, i] <- sqrt(t(part.derivs$Eta)%*%sigma.pars[[i]]%*%part.derivs$Eta)
  }

  ses <- list("se.nie" = se.nie, "se.nde" = se.nde, "se.te" = se.te)
  return(ses)

}

#'@rdname stderrs
#'@export
stderr.cc <- function(Rho, betas, thetas, sigma.pars,
                      x.med, x.out, alt.decomposition, exp.value, control.value){

  se.nie <- matrix(nrow = 1, ncol = length(Rho))
  colnames(se.nie) <- paste(Rho)
  se.nde <- se.nie
  se.te <- se.nie

  t.de <- ifelse(alt.decomposition == TRUE, exp.value, control.value)
  t.ie <- ifelse(alt.decomposition == TRUE, control.value, exp.value)

  for(i in 1:length(Rho)){
    b0 <- betas$beta0[i]
    b1 <- betas$beta1[i]
    b2 <- betas$beta2[, i]
    b3 <- betas$beta3[, i]

    th2 <- thetas$theta2[i]
    th3 <- thetas$theta3[i]
    th6 <- thetas$theta6[, i]
    th7 <- thetas$theta7[, i]

    part.derivs <- partdevs.cc(beta0 = b0, beta1 = b1, beta2 = b2, beta3 = b3, theta2 = th2, theta3 = th3,
                               theta6 = th6, theta7 = th7, exp.value = exp.value, control.value = control.value,
                               x.med = x.med, x.out = x.out, t.de = t.de, t.ie = t.ie)

    se.nie[, i] <- sqrt(t(part.derivs$Gamma)%*%sigma.pars[[i]]%*%part.derivs$Gamma)
    se.nde[, i] <- sqrt(t(part.derivs$Lambda)%*%sigma.pars[[i]]%*%part.derivs$Lambda)
    se.te[, i] <- sqrt(t(part.derivs$Eta)%*%sigma.pars[[i]]%*%part.derivs$Eta)

  }

  ses <- list("se.nie" = se.nie, "se.nde" = se.nde, "se.te" = se.te)
  return(ses)

}


#'Implementations of the partial derivatives (gradients) of the expressions for the direct, indirect and total effects. Used to calculate standard errors (delta method).
#'
#'Functions implementing the partial derivatives (gradients) of the expressions for the direct, indirect and total effects. These are then used
#'to calculate standard errors of the effects using the delta method. Called by the \code{\link{stderrs}} functions. The functions are named according to the convention
#'\code{partdevs."mediator model type""outcome model type"} where \code{b}stands for binary probit regression and
#'\code{c} stands for linear regression.
#'
#'@param beta0,beta1 Vectors of mediator regression parameters (intercept and exposure) over \code{Rho}
#'@param beta2,beta3 Matrices of mediator regression parameters (covariate main effects and exposure-covariate interactions) over \code{Rho}
#'@param theta0,theta1,theta2,theta3 Vectors of outcome regression parameters (intercept, exposure, mediator, exposure-mediator interaction) over \code{Rho}
#'@param theta4,theta5,theta6,theta7 Matrices of outcome regression parameters (covariate main effects, exposure-covariate, mediator-covariate and exposure-mediator-covariate interactions) over \code{Rho}
#'@param sigma.eta For a continuous mediator and binary outcome, matrix with the estimated residual standard deviation for the mediator model over the range of \code{Rho}.
#'@param x.med Mediator covariate matrix for which to calculate standard errors
#'@param x.out Outcome covariate matrix for which to calculate standard errors
#'@param exp.value value of the exposure variable used as the exposure condition.
#'@param control.value value of the exposure variable used as the control (unexposed) condition.
#'@param t.de,t.ie exposure values used to calculate the direct and indirect effects depending on the desired decomposition (see the Details section of \code{\link{sensmediation}} for more information). If \code{alt.decomposition = TRUE} then \code{t.de = exp.value} and \code{t.ie = control.value}, otherwise \code{t.de = control.value} and \code{t.ie = exp.value}.
#'@name partdevs
NULL

#'@rdname partdevs
#'@export
partdevs.bb <-  function(beta0, beta1, beta2, beta3, theta0, theta1, theta2, theta3, theta4, theta5, theta6, theta7,
                         x.med, x.out, t.de, t.ie, exp.value, control.value){

  A.tr <-  theta0 + theta1*exp.value + x.out%*%(theta4 + exp.value*theta5)
  A.cont <-  theta0 + theta1*control.value + x.out%*%(theta4 + control.value*theta5)
  A.ie <-  theta0 + theta1*t.ie + x.out%*%(theta4 + t.ie*theta5)

  B.tr <-  theta0 + theta2 + exp.value*(theta1 + theta3) + x.out%*%(theta4 + exp.value*theta5 +
                                                                      theta6 + exp.value*theta7)
  B.cont <-  theta0 + theta2 + control.value*(theta1 + theta3) + x.out%*%(theta4 + control.value*theta5 +
                                                                            theta6 + control.value*theta7)
  B.ie <-  theta0 + theta2 + t.ie*(theta1 + theta3) + x.out%*%(theta4 + t.ie*theta5 + theta6 + t.ie*theta7)

  C.de <-  beta0 + beta1*t.de + x.med%*%(beta2 + beta3*t.de)
  C.tr <-  beta0 + beta1*exp.value + x.med%*%(beta2 + beta3*exp.value)
  C.cont <-  beta0 + beta1*control.value + x.med%*%(beta2 + beta3*control.value)

  ### NDE ###
  D1 <-  stats::dnorm(C.de)*(stats::pnorm(B.tr) - stats::pnorm(B.cont)) - stats::dnorm(C.de)*(stats::pnorm(A.tr) - stats::pnorm(A.cont))
  d1 <-  mean(D1)

  d2 <-  d1*t.de

  D3 <-  as.vector(D1)*x.med
  d3 <-  colMeans(D3)

  d4 <-  d3*t.de

  D5 <-  (stats::dnorm(A.tr) - stats::dnorm(A.cont))*(1 - stats::pnorm(C.de))  +  (stats::dnorm(B.tr) - stats::dnorm(B.cont))*stats::pnorm(C.de)
  d5 <-  mean(D5)

  D6 <-  (exp.value*stats::dnorm(A.tr) - control.value*stats::dnorm(A.cont))*(1 - stats::pnorm(C.de))  +
    (exp.value*stats::dnorm(B.tr) - control.value*stats::dnorm(B.cont))*stats::pnorm(C.de)
  d6 <-  mean(D6)

  D7 <-  (stats::dnorm(B.tr) - stats::dnorm(B.cont))*stats::pnorm(C.de)
  d7 <-  mean(D7)

  D8 <-  (exp.value*stats::dnorm(B.tr) - control.value*stats::dnorm(B.cont))*stats::pnorm(C.de)
  d8 <-  mean(D8)

  d9 <-  colMeans(as.vector(D5)*x.out)

  d10 <-  colMeans(as.vector(D6)*x.out)

  d11 <-  colMeans(as.vector(D7)*x.out)

  d12 <-  colMeans(as.vector(D8)*x.out)

  Lambda <-  c(d5, d6, d7, d8, d9, d10, d11, d12, d1, d2, d3, d4)

  ### NIE ###
  G1 <-  (stats::pnorm(B.ie) - stats::pnorm(A.ie))*(stats::dnorm(C.tr) - stats::dnorm(C.cont))
  g1 <-  mean(G1)

  G2 <-  (stats::pnorm(B.ie) - stats::pnorm(A.ie))*(exp.value*stats::dnorm(C.tr) - control.value*stats::dnorm(C.cont))
  g2 <-  mean(G2)

  g3 <-  colMeans(as.vector(G1)*x.med)

  g4 <-  colMeans(as.vector(G2)*x.med)

  G5 <-  (stats::dnorm(B.ie) - stats::dnorm(A.ie))*(stats::pnorm(C.tr) - stats::pnorm(C.cont))
  g5 <-  mean(G5)

  g6 <-  t.ie*g5

  G7 <-  stats::dnorm(B.ie)*(stats::pnorm(C.tr) - stats::pnorm(C.cont))
  g7 <-  mean(G7)

  g8 <-  g7*t.ie

  g9 <-  colMeans(as.vector(G5)*x.out)

  g10 <-  t.ie*g9

  g11 <-  colMeans(as.vector(G7)*x.out)

  g12 <-  g11*t.ie

  Gamma <-  c(g5, g6, g7, g8, g9, g10, g11, g12, g1, g2, g3, g4)

  ### TE ###
  out.te.tr1 <-  theta0 + theta1*exp.value + x.out%*%(theta4 + exp.value*theta5)
  out.te.tr2 <-  theta0 + theta2 + exp.value*(theta1 + theta3) + x.out%*%(theta4 + exp.value*theta5 + theta6 +
                                                                            exp.value*theta7)
  out.te.cont1 <-  theta0 + theta1*control.value + x.out%*%(theta4 + control.value*theta5)
  out.te.cont2 <-  theta0 + theta2 + control.value*(theta1 + theta3) + x.out%*%(theta4 + control.value*theta5 +
                                                                                  theta6 + control.value*theta7)
  med.te.tr <-  beta0 + beta1*exp.value + x.med%*%(beta2 + beta3*exp.value)
  med.te.cont <-  beta0 + beta1*control.value + x.med%*%(beta2 + beta3*control.value)

  H1 <-   - stats::pnorm(out.te.tr1)*stats::dnorm(med.te.tr)  +  stats::pnorm(out.te.tr2)*stats::dnorm(med.te.tr)  +
    stats::pnorm(out.te.cont1)*stats::dnorm(med.te.cont)  -  stats::pnorm(out.te.cont2)*stats::dnorm(med.te.cont)
  h1 <-  mean(H1)

  H2 <-  exp.value*( - stats::pnorm(out.te.tr1)*stats::dnorm(med.te.tr)  +  stats::pnorm(out.te.tr2)*stats::dnorm(med.te.tr))  +
    control.value*(stats::pnorm(out.te.cont1)*stats::dnorm(med.te.cont)  -  stats::pnorm(out.te.cont2)*stats::dnorm(med.te.cont))
  h2 <-  mean(H2)

  h3 <-  colMeans(as.vector(H1)*x.med)

  h4 <-  colMeans(as.vector(H2)*x.med)

  H5 <-  stats::dnorm(out.te.tr1)*(1 - stats::pnorm(med.te.tr))  +  stats::dnorm(out.te.tr2)*stats::pnorm(med.te.tr)  -
    stats::dnorm(out.te.cont1)*(1 - stats::pnorm(med.te.cont))  -  stats::dnorm(out.te.cont2)*stats::pnorm(med.te.cont)
  h5 <-  mean(H5)

  H6 <-  exp.value*(stats::dnorm(out.te.tr1)*(1 - stats::pnorm(med.te.tr))  +  stats::dnorm(out.te.tr2)*stats::pnorm(med.te.tr))  -
    control.value*(stats::dnorm(out.te.cont1)*(1 - stats::pnorm(med.te.cont)) - stats::dnorm(out.te.cont2)*stats::pnorm(med.te.cont))
  h6 <-  mean(H6)

  H7 <-  stats::dnorm(out.te.tr2)*stats::pnorm(med.te.tr) - stats::dnorm(out.te.cont2)*stats::pnorm(med.te.cont)
  h7 <-  mean(H7)

  H8 <-  exp.value*stats::dnorm(out.te.tr2)*stats::pnorm(med.te.tr) - control.value*stats::dnorm(out.te.cont2)*stats::pnorm(med.te.cont)
  h8 <-  mean(H8)

  h9 <-  colMeans(as.vector(H5)*x.out)

  h10 <-  colMeans(as.vector(H6)*x.out)

  h11 <-  colMeans(as.vector(H7)*x.out)

  h12 <-  colMeans(as.vector(H8)*x.out)

  Eta <-  c(h5, h6, h7, h8, h9, h10, h11, h12, h1, h2, h3, h4)

  return(list("Lambda" = Lambda, "Gamma" = Gamma, "Eta" = Eta))
}

#'@rdname partdevs
#'@export
partdevs.bc <- function(beta0, beta1, beta2, beta3, theta2, theta3, theta6, theta7, x.med, x.out, t.de, t.ie,
                        exp.value, control.value){

  diff <- exp.value - control.value

  A <- beta0  +  beta1*t.de + x.med%*%(beta2  +  t.de*beta3)
  out.ie <- theta2  +  theta3*t.ie + x.out%*%(theta6  +  t.ie*theta7)
  med.ie <- stats::pnorm(beta0 + beta1*exp.value + x.med%*%(beta2  +  exp.value*beta3))-
    stats::pnorm(beta0 + beta1*control.value + x.med%*%(beta2  +  control.value*beta3))

  ###NDE###
  D1 <- (theta3*diff + x.out%*%theta7*diff)*stats::dnorm(A)
  d1 <- mean(D1)

  d2 <- d1*t.de

  D3 <- as.vector(D1)*x.med
  d3 <- colMeans(D3)

  d4 <- d3*t.de

  d5 <- 0

  d6 <- diff

  d7 <- 0

  D8 <- diff*stats::pnorm(A)
  d8 <- mean(D8)

  d9 <- rep(0, ncol(x.out))

  d10 <- colMeans(diff*x.out)

  d11 <- d9

  d12 <- colMeans(as.vector(D8)*x.out)

  Lambda <- c(d5, d6, d7, d8, d9, d10, d11, d12, d1, d2, d3, d4)

  ###NIE###
  G1 <- out.ie*(stats::dnorm(beta0 + beta1*exp.value + x.med%*%(beta2  +  exp.value*beta3))-
                  stats::dnorm(beta0 + beta1*control.value + x.med%*%(beta2  +  control.value*beta3)))
  g1 <- mean(G1)

  G2 <- out.ie*(exp.value*stats::dnorm(beta0 + beta1*exp.value + x.med%*%(beta2  +  exp.value*beta3))-
                  control.value*stats::dnorm(beta0 + beta1*control.value + x.med%*%(beta2  +  control.value*beta3)))
  g2 <- mean(G2)

  g3 <- colMeans(as.vector(G1)*x.med)

  g4 <- colMeans(as.vector(G2)*x.med)

  g5 <- 0

  g6 <- 0

  G7 <- med.ie
  g7 <- mean(G7)

  g8 <- g7*t.ie

  g9 <- rep(0, ncol(x.out))

  g10 <- g9

  g11 <- colMeans(as.vector(G7)*x.out)

  g12 <- g11*t.ie

  Gamma <- c(g5, g6, g7, g8, g9, g10, g11, g12, g1, g2, g3, g4)

  #TE
  out.te.exp <- theta2  +  theta3*exp.value + x.out%*%(theta6  +  exp.value*theta7)
  out.te.contr <- theta2  +  theta3*control.value + x.out%*%(theta6  +  control.value*theta7)
  med.te.exp <- beta0 + beta1*exp.value + x.med%*%(beta2  +  exp.value*beta3)
  med.te.contr <- beta0 + beta1*control.value + x.med%*%(beta2  +  control.value*beta3)

  H1 <- out.te.exp*stats::dnorm(med.te.exp) - out.te.contr*stats::dnorm(med.te.contr)
  h1 <- mean(H1)

  H2 <- exp.value*out.te.exp*stats::dnorm(med.te.exp) - control.value*out.te.contr*stats::dnorm(med.te.contr)
  h2 <- mean(H2)

  h3 <- colMeans(as.vector(H1)*x.med)

  h4 <- colMeans(as.vector(H2)*x.med)

  h5 <- 0

  h6 <- diff

  H7 <- stats::pnorm(med.te.exp)-stats::pnorm(med.te.contr)
  h7 <- mean(H7)

  H8 <- exp.value*stats::pnorm(med.te.exp)-control.value*stats::pnorm(med.te.contr)
  h8 <- mean(H8)

  h9 <- g9

  h10 <- colMeans(diff*x.out)

  h11 <- colMeans(as.vector(H7)*x.out)

  h12 <- colMeans(as.vector(H8)*x.out)

  Eta <- c(h5, h6, h7, h8, h9, h10, h11, h12, h1, h2, h3, h4)

  return(list("Lambda"=Lambda, "Gamma"=Gamma, "Eta"=Eta))
}

#'@rdname partdevs
#'@export
partdevs.cb <- function(beta0, beta1, beta2, beta3, theta0, theta1, theta2, theta3, theta4,
                        theta5, theta6, theta7, sigma.eta, x.med, x.out, t.de, t.ie, exp.value, control.value){

  med.ie.tr <- beta0 + beta1*exp.value + x.med%*%(beta2 + exp.value*beta3)
  med.ie.cont <- beta0 + beta1*control.value + x.med%*%(beta2 + control.value*beta3)
  med.de <- beta0 + beta1*t.de + x.med%*%(beta2 + t.de*beta3)

  out.ie <- theta2 + theta3*t.ie + x.out%*%(theta6 + t.ie*theta7)
  out.de.tr <- theta2 + theta3*exp.value + x.out%*%(theta6 + exp.value*theta7)
  out.de.cont <- theta2 + theta3*control.value + x.out%*%(theta6 + control.value*theta7)

  denom.ie <- sqrt(sigma.eta^2*out.ie^2 + 1)
  denom.de.tr <- sqrt(sigma.eta^2*out.de.tr^2 + 1)
  denom.de.cont <- sqrt(sigma.eta^2*out.de.cont^2 + 1)

  C.tr <- theta0 + theta1*exp.value + x.out%*%(theta4 + theta5*exp.value)
  C.cont <- theta0 + theta1*control.value + x.out%*%(theta4 + theta5*control.value)
  C.ie <- theta0 + theta1*t.ie + x.out%*%(theta4 + theta5*t.ie)

  ###NDE###
  D1 <- out.de.tr/denom.de.tr*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    out.de.cont/denom.de.cont*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d1 <- mean(D1)

  d2 <- d1*t.de

  D3 <- as.vector(D1)*x.med
  d3 <- colMeans(D3)

  d4 <- d3*t.de

  D5 <- -(C.tr + out.de.tr*med.de)*sigma.eta*out.de.tr^2/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr ) +
    (C.cont + out.de.cont*med.de)*sigma.eta*out.de.cont^2/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d5 <- mean(D5)

  D6 <- 1/denom.de.tr*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    1/denom.de.cont*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d6 <- mean(D6)

  D7 <- 1/denom.de.tr*exp.value*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    1/denom.de.cont*control.value*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d7 <- mean(D7)

  D8 <- (med.de - C.tr*sigma.eta^2*out.de.tr)/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    (med.de - C.cont*sigma.eta^2*out.de.cont)/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d8 <- mean(D8)

  D9 <-  exp.value*(med.de - C.tr*sigma.eta^2*out.de.tr)/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.de +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    control.value*(med.de - C.cont*sigma.eta^2*out.de.cont)/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.de +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  d9 <- mean(D9)

  d10 <- colMeans(as.vector(D6)*x.out)

  d11 <- colMeans(as.vector(D7)*x.out)

  d12 <- colMeans(as.vector(D8)*x.out)

  d13 <- colMeans(as.vector(D9)*x.out)

  Lambda <- c(d6, d7, d8, d9, d10, d11, d12, d13, d1, d2, d3, d4, d5)

  ###NIE###
  G1 <- out.ie/denom.ie*(stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.tr +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )  -
                           stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.cont +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie ))
  g1 <- mean(G1)

  G2 <- out.ie/denom.ie*(exp.value*stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.tr +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )  -
                           control.value*stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.cont +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie ))
  g2 <- mean(G2)

  g3 <- colMeans(as.vector(G1)*x.med)

  g4 <- colMeans(as.vector(G2)*x.med)

  G5 <- sigma.eta*out.ie^2/(denom.ie^3)*(-stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.tr +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )*(C.ie + out.ie*med.ie.tr) +
                                           stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.cont +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie ) *(C.ie + out.ie*med.ie.cont) )
  g5 <- mean(G5)

  G6 <- 1/denom.ie*(stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.tr +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )  -
                      stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.cont +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie ))
  g6 <- mean(G6)

  g7 <- t.ie*g6

  G8 <- (med.ie.tr - C.ie*sigma.eta^2*out.ie)/(denom.ie^3)*stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.tr +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )  -
    (med.ie.cont - C.ie*sigma.eta^2*out.ie)/(denom.ie^3)*stats::dnorm((theta0 + theta1*t.ie + out.ie*med.ie.cont +  x.out%*%(theta4 + theta5*t.ie) )/denom.ie )
  g8 <- mean(G8)

  g9 <- g8*t.ie

  g10 <- colMeans(as.vector(G6)*x.out)

  g11 <- colMeans(as.vector(G6*t.ie)*x.out)

  g12 <- colMeans(as.vector(G8)*x.out)

  g13 <- colMeans(as.vector(G8*t.ie)*x.out)

  Gamma <- c(g6, g7, g8, g9, g10, g11, g12, g13, g1, g2, g3, g4, g5)

  ###TE###
  H1 <- out.de.tr/denom.de.tr*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    out.de.cont/denom.de.cont*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h1 <- mean(H1)

  H2 <- exp.value*out.de.tr/denom.de.tr*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    control.value*out.de.cont/denom.de.cont*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h2 <- mean(H2)

  h3 <- colMeans(as.vector(H1)*x.med)

  h4 <- colMeans(as.vector(H2)*x.med)

  H5 <- -(C.tr + out.de.tr*med.ie.tr)*sigma.eta*out.de.tr^2/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr ) +
    (C.cont + out.de.cont*med.ie.cont)*sigma.eta*out.de.cont^2/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h5 <- mean(H5)

  H6 <- 1/denom.de.tr*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    1/denom.de.cont*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h6 <- mean(H6)

  H7 <- 1/denom.de.tr*exp.value*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    1/denom.de.cont*control.value*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h7 <- mean(H7)

  H8 <- (med.ie.tr - C.tr*sigma.eta^2*out.de.tr)/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    (med.ie.cont - C.cont*sigma.eta^2*out.de.cont)/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h8 <- mean(H8)

  H9 <-  exp.value*(med.ie.tr - C.tr*sigma.eta^2*out.de.tr)/(denom.de.tr^3)*stats::dnorm((theta0 + theta1*exp.value + out.de.tr*med.ie.tr +  x.out%*%(theta4 + theta5*exp.value) )/denom.de.tr )  -
    control.value*(med.ie.cont - C.cont*sigma.eta^2*out.de.cont)/(denom.de.cont^3)*stats::dnorm((theta0 + theta1*control.value + out.de.cont*med.ie.cont +  x.out%*%(theta4 + theta5*control.value) )/denom.de.cont )
  h9 <- mean(H9)

  h10 <- colMeans(as.vector(H6)*x.out)

  h11 <- colMeans(as.vector(H7)*x.out)

  h12 <- colMeans(as.vector(H8)*x.out)

  h13 <- colMeans(as.vector(H9)*x.out)

  Eta <- c(h6, h7, h8, h9, h10, h11, h12, h13, h1, h2, h3, h4, h5)

  return(list("Lambda"=Lambda, "Gamma"=Gamma, "Eta"=Eta))
}

#'@rdname partdevs
#'@export
partdevs.cc <- function(beta0, beta1, beta2, beta3, theta2, theta3, theta6, theta7, exp.value, control.value,
                        x.med, x.out, t.de, t.ie){

  diff <- exp.value - control.value
  A <- beta0 + beta1*t.de + x.med%*%(beta2 + t.de*beta3)
  B <- theta2 + theta3*t.ie + x.out%*%(theta6 + t.ie*theta7)

  ###NDE###
  D1 <- theta3*diff + x.out%*%theta7*diff
  d1 <- mean(D1)

  d2 <- d1*t.de

  D3 <- as.vector(D1)*x.med
  d3 <- colMeans(D3)

  d4 <- d3*t.de

  d5 <- 0

  d6 <- diff

  d7 <- 0

  D8 <- diff*A
  d8 <- mean(D8)

  d9 <- rep(0, ncol(x.out))

  d10 <- colMeans(diff*x.out)

  d11 <- d9

  d12 <- colMeans(as.vector(D8)*x.out)

  Lambda <- c(d5, d6, d7, d8, d9, d10, d11, d12, d1, d2, d3, d4)

  ###NIE###
  g1 <- 0

  G2 <- diff*B
  g2 <- mean(G2)

  g3 <- rep(0, ncol(x.med))

  g4 <- colMeans(as.vector(G2)*x.med)

  g5 <- 0

  g6 <- 0

  G7 <- diff*(beta1 + x.med%*%beta3)
  g7 <- mean(G7)

  g8 <- g7*t.ie

  g9 <- rep(0, ncol(x.out))

  g10 <- g9

  g11 <- colMeans(as.vector(G7)*x.out)

  g12 <- g11*t.ie

  Gamma <- c(g5, g6, g7, g8, g9, g10, g11, g12, g1, g2, g3, g4)

  ###TE###
  med.te <- (beta1 + x.med%*%beta3)*diff
  med.te2 <- exp.value*(beta0 + beta1*exp.value + x.med%*%(beta2 + exp.value*beta3))-
    control.value*(beta0 + beta1*control.value + x.med%*%(beta2 + control.value*beta3))

  H1 <- theta3*diff + x.out%*%theta7*diff
  h1 <- mean(H1)

  H2 <- theta2*diff + theta3*(exp.value^2-control.value^2) + x.out%*%theta6*diff +
    x.out%*%theta7*(exp.value^2-control.value^2)
  h2 <- mean(H2)

  h3 <- colMeans(as.vector(H1)*x.med)

  h4 <- colMeans(as.vector(H2)*x.med)

  h5 <- 0

  h6 <- diff

  H7 <- med.te
  h7 <- mean(H7)

  H8 <- med.te2
  h8 <- mean(H8)

  h9 <- g9

  h10 <- colMeans(h6*x.out)

  h11 <- colMeans(as.vector(H7)*x.out)

  h12 <- colMeans(as.vector(H8)*x.out)

  Eta <- c(h5, h6, h7, h8, h9, h10, h11, h12, h1, h2, h3, h4)

  return(list("Lambda"=Lambda, "Gamma"=Gamma,  "Eta"=Eta))
}

