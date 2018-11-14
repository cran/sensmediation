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
