# sensmediation 0.3.0

* Fixed a bug in the `calc.effects` function that could affect the estimated direct and indirect effects with moderation on more than one covariate
if the exposure-covariate, mediator-covariate, exposure-mediator-covariate interaction terms were not in the same order as the main effects of the covariates
* Fixed an error in the log-likelihoods for cases with probit `model.expl` and continuous `model.resp` as well as cases with probit `model.resp` 
and continuous `model.expl` which gave results for -rho instead of rho and vice versa.
* Corrected warning messages for omitted `exp.name` and `med.name` arguments in `sensmediation`
* Estimated outcome and mediator model parameters separated into different categories depending on effect type are now stored in `effectsMed` objects
* The partial derivatives of the NDE, NIE and TE wrt the mediator and outcome model parameters for each value of Rho are now stored in `effectsMed` objects
* A list with the joint covariance matrix of the outcome and mediator model parameters for each value of Rho is now stored in `effectsMed` objects
* A data set consisting of a subsample from Riksstroke, the Swedish Stroke Register, is now provided and used to exemplify the package functions.
* The covariate names in the list of the `covariates` argument now can (need to) be given as coded in the original data rather than the way they are recoded in 
the glm output.
* A note is printed if one or more covariates in `covariates` are not found in the mediator and outcome models and therefore not conditioned on. Variables not conditioned on are not printed
in the output from summary.effectsMed.

# sensmediation 0.2.0

* Optimization function changed from `optimx` to `maxLik`
* Analytic Hessians of the log-likelihoods have been implemented
* Default optimization method changed from BFGS to Newton Raphson
* Arguments allowed to be passed to the maximization function (`method` and `control`)
* The `sensmediation` argument `covariance` has been deprecated
* The `sensmediation` arguments `out.full`, `med.full` and `all.interactions` have been deprecated
* New `sensmediation` arguments: `exp.name` and `med.name`
* Errors in the analytic gradients of the cases with probit `model.expl` and continuous `model.resp` as well as cases with probit `model.resp` 
and continuous `model.expl` have been fixed.
