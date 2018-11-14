#'Summary function for objects of class \code{"effectsMed"}
#'
#'@param object object of class \code{"effectsMed"}
#'@param non.sign logical indicating whether sensitivity analysis results should be printed for non-significant effects.
#'@param ... additional arguments
#'@return A list with values:
#'\item{call}{The matched call}
#'\item{Rho}{The sensitivity parameter vector.}
#'\item{type}{character, the type of confounding the sensitivity analysis is performed for.}
#'\item{conf.level}{numeric, the confidence level used for confidence intervals and uncertainty intervals.}
#'\item{UI}{matrix with the uncertainty intervals for the NIE (NIE*) and NDE (NDE*) over the range of the sensitivity parameter \code{Rho}.}
#'\item{covariates}{list of the covariate values that the effects are conditioned on.}
#'\item{exp.name}{character vector containing the name of the exposure.}
#'\item{med.name}{character vector containing the name of the mediator.}
#'\item{alt.decomposition}{logical, indicating whether the alternative definitions of the direct and indirect effects have been used}
#'\item{non.sign}{logical indicating whether sensitivity analysis results are printed for non-significant effects.}
#'\item{effects}{Results of the mediation analysis. Estimated NIE and NDE with confidence intervals and p-values for \code{Rho = 0}}
#'\item{ns.nie}{values of \code{Rho} with estimated NIEs and confidence intervals where the NIE is not significant.}
#'\item{ns.nde}{values of \code{Rho} with estimated NDEs and confidence intervals where the NDE is not significant.}
#'\item{rev.nie}{values of \code{Rho} with estimated NIEs and confidence intervals where the NIE is reversed.}
#'\item{rev.nde}{values of \code{Rho} with estimated NDEs and confidence intervals where the NDE is reversed.}
#'@export

summary.effectsMed <- function(object, non.sign = FALSE,...)
{
	x <- object

	ans <- list()
	ans$call <- x$call
	ans$Rho <- x$Rho
	ans$conf.level <- x$conf.level
	ans$UI <- x$UI
	ans$type <- x$type
	ans$covariates <- x$covariates
	ans$exp.name <- x$exp.name
	ans$med.name <- x$med.name
	ans$non.sign <- non.sign
	ans$alt.decomposition <- x$alt.decomposition


	i0 <- which(x$Rho == 0)
	est <- as.matrix(c(x$NIE[i0], x$NDE[i0], x$NIE[i0] + x$NDE[i0]))
	se <- as.matrix(c(x$std.errs$se.nie[i0], x$std.errs$se.nde[i0], x$std.errs$se.te[i0]))
	zval <- est/se
	lower <- as.matrix(c(x$CI$CI.nie[i0, 1], x$CI$CI.nde[i0, 1], x$CI$CI.te[i0, 1]))
	upper <- as.matrix(c(x$CI$CI.nie[i0, 2], x$CI$CI.nde[i0, 2], x$CI$CI.te[i0, 2]))

	ans$effects <- cbind(est, lower, upper, 2 * stats::pnorm(abs(zval), lower.tail = FALSE))
	colnames(ans$effects) <- c("Estimate", gettextf("%d%% CI Lower",(100*x$conf.level)), gettextf("%d%% CI Upper",(100*x$conf.level)),
	                           "p-value")
	if(x$alt.decomposition == TRUE)
	  rownames(ans$effects) <- c("NIE*","NDE*","TE")
	else
	  rownames(ans$effects) <- c("NIE","NDE","TE")

	pos.NS.nie <- which(sign(x$CI$CI.nie[,1]) != sign(x$CI$CI.nie[,2]))
	ans$NS.nie <- cbind(x$NIE[pos.NS.nie], x$CI$CI.nie[pos.NS.nie, 1], x$CI$CI.nie[pos.NS.nie, 2])
	colnames(ans$NS.nie) <- c("Estimate", gettextf("%d%% CI Lower",(100*x$conf.level)), gettextf("%d%% CI Upper",(100*x$conf.level)))
	rownames(ans$NS.nie) <- paste(x$Rho[pos.NS.nie])

	pos.NS.nde <- which(sign(x$CI$CI.nde[,1])!=sign(x$CI$CI.nde[,2]))
	ans$NS.nde <- cbind(x$NDE[pos.NS.nde], x$CI$CI.nde[pos.NS.nde, 1], x$CI$CI.nde[pos.NS.nde, 2])
	colnames(ans$NS.nde) <- c("Estimate", gettextf("%d%% CI Lower",(100*x$conf.level)), gettextf("%d%% CI Upper",(100*x$conf.level)))
	rownames(ans$NS.nde) <- paste(x$Rho[pos.NS.nde])


	pos.rev.nie <- which(sign(x$CI$CI.nie[,1])==-sign(est[1,])&sign(x$CI$CI.nie[,2])==-sign(est[1,]))
	ans$rev.nie <- cbind(x$NIE[pos.rev.nie], x$CI$CI.nie[pos.rev.nie, 1], x$CI$CI.nie[pos.rev.nie, 2])
	colnames(ans$rev.nie) <- c("Estimate", gettextf("%d%% CI Lower",(100*x$conf.level)), gettextf("%d%% CI Upper",(100*x$conf.level)))
	rownames(ans$rev.nie) <- paste(x$Rho[pos.rev.nie])

	pos.rev.nde <- which(sign(x$CI$CI.nde[,1]) == -sign(est[2,]) & sign(x$CI$CI.nde[,2]) == -sign(est[2,]))
	ans$rev.nde <- cbind(x$NDE[pos.rev.nde], x$CI$CI.nde[pos.rev.nde, 1], x$CI$CI.nde[pos.rev.nde, 2])
	colnames(ans$rev.nde) <- c("Estimate", gettextf("%d%% CI Lower",(100*x$conf.level)), gettextf("%d%% CI Upper",(100*x$conf.level)))
	rownames(ans$rev.nde) <- paste(x$Rho[pos.rev.nde])

	class(ans) <- c("summaryeffectsMed", "effectsMed")
	ans

}

#'Plot function for objects of class \code{"effectsMed"}
#'
#'Plots the estimated natural indirect or direct effects with confidence intervals over the range of the sensitivity parameter \code{Rho}.
#'@param x object of class \code{"effectsMed"}
#'@param effect which effect to plot results for ("indirect" or "direct")
#'@param xlab a title for the x axis, see \code{\link{title}}. Default is \code{expression(rho)}.
#'@param ylab a title for the y axis, see \code{\link{title}}. Default is \code{NIE} (\code{NIE*} if \code{alt.decomposition = TRUE}) or \code{NDE} (\code{NDE*})
#'@param xlim the x limits (x1, x2) of the plot, see \code{\link{plot.default}}. Default is \code{c(min(x$Rho),max(x$Rho))}
#'@param ylim the y limits of the plot. Default is \code{c(min(x$CI$CI.nie[,1]),max(x$CI$CI.nie[,2]))}
#'@param main a main title for the plot, see \code{\link{title}}
#'@param lwd line widths for the lines of the plot, see \code{\link{par}}
#'@param ... additional graphical parameters to be passed to plotting functions, see \code{\link{par}}
#'@export
#'

plot.effectsMed <- function(x, effect="indirect", xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL, lwd = graphics::par("lwd"),...){

  if(is.null(xlab))
    xlab <- expression(rho)
  if(is.null(ylab)){
    if(effect == "indirect"){
      if(x$alt.decomposition == TRUE)
        ylab <- "NIE*"
      else
        ylab <- "NIE"
    }
    else{
      if(x$alt.decomposition == TRUE)
        ylab <- "NDE*"
      else
        ylab <- "NDE"
    }

  }


  if(effect == "indirect"){
    lower <- x$CI$CI.nie[,1]
    upper <- x$CI$CI.nie[,2]

    if(is.null(ylim))
      ylim <- c(min(lower), max(upper))

    if(is.null(xlim))
      xlim <- c(min(x$Rho), max(x$Rho))

    graphics::plot(x$Rho, x$NIE, xlab = xlab, ylab = ylab, type = "l", ylim = ylim, xlim = xlim,...)
    graphics::polygon(c(x$Rho, rev(x$Rho)), c(lower, rev(upper)), col = "grey80", border = FALSE)
    graphics::lines(x$Rho, x$NIE, lwd = lwd)
    graphics::abline(h = 0, lty = 5, lwd = lwd)
    graphics::points(0, x$NIE[which(x$Rho == 0)], pch = 20)
    ranges <- graphics::par("usr")
    graphics::segments(ranges[1],x$NIE[which(x$Rho == 0)], 0, x$NIE[which(x$Rho == 0)], lty = 3, lwd = lwd)
    graphics::segments(0, x$NIE[which(x$Rho == 0)], 0, ranges[3], lty = 3, lwd = lwd)

  }

  if(effect == "direct"){
    lower <- x$CI$CI.nde[,1]
    upper <- x$CI$CI.nde[,2]

    if(is.null(ylim))
      ylim <- c(min(lower), max(upper))

    if(is.null(xlim))
      xlim <- c(min(x$Rho), max(x$Rho))


    graphics::plot(x$Rho, x$NDE, xlab = xlab, ylab = ylab, type = "l", ylim = ylim, xlim = xlim,...)
    graphics::polygon(c(x$Rho, rev(x$Rho)), c(lower, rev(upper)), col = "grey80", border = FALSE)
    graphics::lines(x$Rho, x$NDE, lwd = lwd)
    graphics::abline(h = 0, lty = 5, lwd = lwd)
    graphics::points(0, x$NDE[which(x$Rho == 0)], pch = 20)
    ranges <- graphics::par("usr")
    graphics::segments(ranges[1],x$NDE[which(x$Rho == 0)], 0, x$NDE[which(x$Rho == 0)],lty = 3, lwd = lwd)
    graphics::segments(0, x$NDE[which(x$Rho == 0)], 0, ranges[3],lty = 3, lwd = lwd)

  }

  if(is.null(main)){
    if(x$type == "my")
      graphics::title(main="Sensitivity analysis mediator-outcome confounding")
    if(x$type == "zm")
      graphics::title(main="Sensitivity analysis exposure-mediator confounding")
    if(x$type == "zy")
      graphics::title(main="Sensitivity analysis exposure-outcome confounding")

  }
  else
    graphics::title(main)
  }

#'Print function for objects of class \code{"summaryeffectsMed"}
#'
#'@param x object of class \code{"summaryeffectsMed"}
#'@param digits number of digits to be printed.
#'@param ... additional arguments
#'@export
#'


print.summaryeffectsMed <- function(x, digits= max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n\n", sep = "")

  cat("Mediation analysis:\n\n")

  cat("Exposure =", x$exp.name, "  Mediator =", x$med.name, "\n\n")

  if(length(x$covariates)){
    cat("Effects conditional on the covariate values\n\n")
    print.default(format(unlist(x$covariates),digits = digits), print.gap = 2.5, quote = FALSE)
    cat("\n\n")
  }
  else
    cat("Marginal effects\n\n")

  if (length(x$effects))
    stats::printCoefmat(x$effects, digits = digits, print.gap = 3)

  else cat("No coefficients\n")

  cat("\n\n\n")

  if((x$non.sign == TRUE | x$effects[1, 4] < 1 - x$conf.level | x$effects[2, 4] < 1 - x$conf.level) & length(which(x$Rho != 0))){
    if(x$type == "zm")
      cat("Sensitivity analysis for exposure-mediator confounding:\n\n")

    if(x$type == "my")
      cat("Sensitivity analysis for mediator-outcome confounding:\n\n")

    if(x$type == "zy")
      cat("Sensitivity analysis for exposure-outcome confounding:\n\n")

    cat(gettextf("At least %d%% uncertainty intervals over Rho from %g to %g:\n", (100*x$conf.level), min(x$Rho), max(x$Rho)))
    print.default(x$UI, digits = digits, quote = FALSE, print.gap = 2)


    cat("\n")
  }

  if((x$non.sign == TRUE | x$effects[1, 4] < 1 - x$conf.level) & length(which(x$Rho != 0))){
    # Sensitivity regions
    if(length(x$NS.nie)){
      if(x$alt.decomposition == TRUE)
        cat("NIE* non-significant for:\n")
      else
        cat("NIE non-significant for:\n")
      print.default(x$NS.nie, digits = digits, quote = FALSE, print.gap = 2)
      cat("\n")
    }

    if(length(x$rev.nie)){
      if(x$alt.decomposition == TRUE)
        cat("NIE* reversed for:\n")
      else
        cat("NIE reversed for:\n")
      print.default(x$rev.nie, digits = digits, quote = FALSE, print.gap = 2)
      cat("\n")
    }
  }


  if((x$non.sign == TRUE | x$effects[2, 4] < 1 - x$conf.level) & length(which(x$Rho != 0))){
    if(length(x$NS.nde)){
      if(x$alt.decomposition == TRUE)
        cat("NDE* non-significant for:\n")
      else
        cat("NDE non-significant for:\n")
      print.default(x$NS.nde, digits = digits, quote = FALSE, print.gap = 2)
      cat("\n")
    }

    if(length(x$rev.nde)){
      if(x$alt.decomposition == TRUE)
        cat("NDE* reversed for:\n")
      else
        cat("NDE reversed for:\n")
      print.default(x$rev.nde, digits = digits, quote = FALSE, print.gap = 2)
      cat("\n")
    }
  }


  invisible(x)

}












