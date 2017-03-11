# Implementation of the POUMM likelihood and heritability estimators

#'@title The Phylogenetic (Ornstein-Uhlenbeck) Mixed Model
#'  
#'@description This is the high-level entry point to the POUMM method. The POUMM
#'  function fits the POUMM method to a tree and observed trait-values at its tips
#'   and returns an object of class "POUMM".
#'  
#'@param z Either a numeric vector containing the phenotypic values at the tips 
#'  of tree or a named list containing named elements z - a numeric vector and 
#'  tree - a phylo object (it is possible to specify different element names 
#'  using the arguments zName and treeName).
#'@param tree A phylo object or NULL in case z is a list.
#'@param zName,treeName Character strings used when the parameter z is a list; 
#'  indicate the names in the list of the values-vector and the tree. Default: 
#'  'z' and 'tree'.
#'@param parDigits Integer specifying rounding to be done on the parameter 
#'  vector before likelihood calculation. Defaults to 6 decimal digits. This can
#'  be useful during maximum likelihood optimization to prevent likelihood 
#'  calculation on very small but positive values of alpha, but should be used 
#'  with caution since specifying a small number of digits, i.e. 2 or 3 can
#'  result in an infinite loop during optim. Specify a negative number
#'  to disable rounding.
#'@param ... additional arguments passed to the likPOUMMGivenTreeVTips function 
#'  (?dVGivenTreeOU for details).
#'@param spec A named list specifying how the ML and MCMC fit should be done. 
#'  See ?specifyPOUMM.
#'@param doMCMC logical: should a MCMC fit be performed. An MCMC fit provides a 
#'  sample from the posterior distribution of the parameters given a prior 
#'  distribution and the data. Unlike the ML-fit, it allows to estimate 
#'  confidence intervals for the estimated parameters. This argument is TRUE by 
#'  default. The current implementation uses a modified version of the adaptive 
#'  Metropolis sampler from the package "adaptMCMC" written by Andreas 
#'  Scheidegger. To obtain meaningful estimates MCMC may need to run for several
#'  millions of iterations (parameter nSamplesMCMC set to 1e5 by default). See 
#'  parameters ending at MCMC in ?specifyPOUMM for details.
#'  
#'@param usempfr integer indicating if and how mpfr should be used for small 
#'  parameter values (any(c(alpha, sigma, sigmae) < 0.01)). Using the mpfr 
#'  package can be forced by specifying an integer greater or equal to 2. 
#'  Setting usempfr=0 disables high precision likelihood calculation. Requires 
#'  the Rmpfr package. Note that when using mpfr, the time for one likelihood 
#'  calculation can increase more than 100-fold. Default (0). Note that during 
#'  the ML and MCMC fit this flag is temporarily raised by 1 in order to enable 
#'  mpfr-check on parameters resulting in higher likelihood than the current 
#'  maximum one.
#'  
#'@param useArma Logical indicating if the Armadillo library should be used for 
#'  faster vector operations. Defaults to TRUE. Since Armadillo doesn't support 
#'  mpfr, it gets disabled when usempfr is bigger than 0.
#'@param verbose,debug Logical flags indicating whether to print informative 
#'  and/or debug information on the standard output (both are set to to FALSE by
#'  default).
#'  
#'@return An object of S3 class 'POUMM'. This object can be analyzed using 
#'  S3 generic functions: \code{\link{summary}}, 
#'  \code{\link{plot}}, \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{coef}},
#'  \code{\link{logLik}}, \code{\link{fitted}}.
#' 
#' @seealso \code{\link{specifyPOUMM}} for parametrizations and custom settings
#'  of the POUMM fit.
#'  
#' @examples 
#' \dontrun{
#' # Please, read the package vignette for more detailed examples.
#' N <- 500
#' tr <- ape::rtree(N)
#' z <- POUMM::rVNodesGivenTreePOUMM(tr, 0, 2, 3, 1, 1)[1:N]
#' fit <- POUMM::POUMM(z, tr, spec = POUMM::specifyPOUMM(nSamplesMCMC = 5e4))
#' plot(fit)
#' summary(fit)
#' AIC(fit)
#' BIC(fit)
#' coef(fit)
#' logLik(fit)
#' fitted(fit)
#' plot(resid(fit))
#' abline(h=0)
#' 
#' # fit PMM to the same data and do a likelihood ratio test
#' fitPMM <- POUMM::POUMM(z, tr, spec = POUMM::specifyPMM(nSamplesMCMC = 5e4))
#' lmtest::lrtest(fitPMM, fit)
#' }
#' 
#'@references 
#'  Mitov, V., and Stadler, T. (2017). POUMM: An R-package for Bayesian Inference 
#'  of Phylogenetic Heritability. bioRxiv, 115089. 
#'  https://doi.org/10.1101/115089
#'  
#'  Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced 
#'  acceptance rate. Statistics and Computing, 22(5), 997-1008. 
#'  http://doi.org/10.1007/s11222-011-9269-5 
#'  
#'  Scheidegger, A. (2012). adaptMCMC: Implementation of a generic adaptive 
#'  Monte Carlo Markov Chain sampler. 
#'  http://CRAN.R-project.org/package=adaptMCMC
#'
#'@importFrom stats var sd rnorm dnorm dexp rexp dunif runif nobs
#'@useDynLib POUMM
#' 
#' @export
POUMM <- function(
  z, tree, zName = 'z', treeName = 'tree', 
  parDigits = 6, usempfr = 0, useArma = TRUE,
  ..., 
  spec = NULL, doMCMC = TRUE,
  verbose = FALSE, debug=FALSE) {
  
  ###### Verify input data ######
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    
    if(is.null(z) | is.null(tree)) {
      stop('If a list is supplied as argument z, this list should contain a ',
           'vector of trait values named "z" or zName and a phylo-object named ',
           '"tree" or treeName')
    }
    if(any(is.na(z)) | any(is.infinite(z))) {
      stop('Check z for infinite or NA values!')
    }
    if(any(tree$edge.length <= 0) | any(is.infinite(tree$edge.length)) | 
       any(is.na(tree$edge.length))) {
      stop('Check the tree for non-finite or non-positive edge-lengths!')
    }
  }
  
  ######## Caching pruneInfo for faster likelihood calculations
  if(!validateZTree(z, tree)) {
    stop("Invalid z and/or tree.")
  }
  
  pruneInfo <- pruneTree(tree, z)
  
  tTips <- nodeTimes(tree, tipsOnly = TRUE)
  
  
  ######## Default POUMM spec ###########
  if(is.function(spec)) {
    spec <- do.call(spec, 
                    list(z = z, tree = tree, 
                         zMin = min(z), zMean = mean(z), zMax = max(z), 
                         zVar = var(z), zSD = sd(z), 
                         tMin = min(tTips), tMean = mean(tTips), tMax = max(tTips)))
  } else {
    spec <- 
      do.call(specifyPOUMM, 
              c(list(z = z, tree = tree, 
                     zMin = min(z), zMean = mean(z), zMax = max(z), 
                     zVar = var(z), tMin = min(tTips), tMean = mean(tTips),
                     tMax = max(tTips)), spec))
  }
  
  result <- list(pruneInfo = pruneInfo, 
                 N = length(tree$tip.label), 
                 tMax = max(nodeTimes(tree, tipsOnly = TRUE)),
                 tMean = mean(nodeTimes(tree, tipsOnly = TRUE)), 
                 spec = spec, 
                 ...)
  
  
  
  # define a loglik function to be called during likelihood maximization as well
  # as MCMC-sampling. 
  # The argument par is a named numeric vector. This vector is mapped to the 
  # POUMM parameters alpha, theta, sigma, sigmae and g0 using the function 
  # spec$parMapping. If g0 = NA the loglik
  # funciton finds the value of g0 that maximizes the likelihood, i.e. 
  # p(z | tree, alpha, theta, sigma, sigmae, g0), or if g0Prior specifiies a 
  # normal distribution with mean g0Prior$mean and variance g0Prior$var, 
  # p(z | tree, alpha, theta, sigma, sigmae, g0) x N(g0 | g0Prior$mean,
  # g0Prior$var). 
  loglik <- function(par, pruneInfo, useCpp = useArma, memo = NULL) {
    if(parDigits >= 0) {
      par <- as.vector(round(par, parDigits))
    }
    
    atsseg0 <- spec$parMapping(par)
    
    if(useCpp & usempfr == 0) {
      # no use of high-precision floating point ops, so we call the faster C++ 
      # likelihood implementation
      val <- likPOUMMGivenTreeVTipsC(
        pruneInfo$integrator, 
        alpha = atsseg0[1], theta = atsseg0[2],
        sigma = atsseg0[3], sigmae = atsseg0[4], g0 = atsseg0[5], 
        g0Prior = spec$g0Prior, log = TRUE)
    } else {
      val <- 
        do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(atsseg0), list(
          g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
          usempfr = usempfr, ...)))
    }
    
    if(is.na(val) | is.infinite(val)) {
      val <- -1e100
      attr(val, "g0") <- atsseg0[5]
      attr(val, "g0LogPrior") <- NA
    } 
    
    if(!is.null(memo)) {
      valMemo <- mget('val', memo, ifnotfound = list(-Inf))$val  
    } else {
      valMemo <- NA
    }
    
    if(!is.na(valMemo) & valMemo + 1 < val) {
      if(usempfr == 0) {
        if(verbose) {
          print(par)
          cat('Rmpfr check on: ', val, '>', valMemo)
        }
        
        val <- 
          do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(atsseg0), list(
            g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
            usempfr = usempfr + 1, ...)))
        
        if(is.na(val) | is.infinite(val)) {
          val <- -1e100
          attr(val, "g0") <- atsseg0[5]
          attr(val, "g0LogPrior") <- NA
        } 
        
        if(verbose) {
          cat('. New: ', val, '.\n')
        }  
      }
    }
    
    val
  }
  
  result$loglik <- loglik
  
  if(verbose) {
    print('Performing ML-fit...')
  }
  fitML <- do.call(
    maxLikPOUMMGivenTreeVTips, 
    c(list(loglik = loglik, verbose = verbose, debug = debug, 
           pruneInfo = pruneInfo), spec))
  
  if(verbose) {
    cat("max loglik from ML: \n")
    print(fitML$value)
    cat("parameters at max loglik from ML: \n")
    print(fitML$par)
  }
  
  result[['fitML']] <- fitML
  
  if(doMCMC) {
    if(verbose) {
      print('Performing MCMC-fit...')
    }
    
    fitMCMC <- do.call(
      mcmcPOUMMGivenPriorTreeVTips, 
      c(list(loglik = loglik, fitML = fitML, verbose = verbose, debug = debug, 
             pruneInfo = pruneInfo), spec))
    
    if(verbose) {
      cat("max loglik from MCMCs: \n")
      print(fitMCMC$valueMaxLoglik)
      cat("parameters at max loglik from MCMCs: \n")
      print(fitMCMC$parMaxLoglik)
    }
    
    # if the max likelihood from the MCMC is bigger than the max likelihood 
    # from the ML-fit, it is likely that the ML fit got stuck in a local 
    # optimum, so we correct it.
    if(fitML$value < fitMCMC$valueMaxLoglik) {
      parInitML <- as.vector(fitMCMC$parMaxLoglik)
      names(parInitML) <- names(spec$parLower)
      
      if(all(c(parInitML >= spec$parLower, parInitML <= spec$parUpper))) {
        if(verbose) {
          cat("The MCMC-fit found a better likelihood than the ML-fit. Performing ML-fit starting from the MCMC optimum.")
        }
        
        spec[["parInitML"]] <- parInitML
        fitML2 <- do.call(
          maxLikPOUMMGivenTreeVTips, 
          c(list(loglik = loglik, verbose = verbose,
                 debug = debug, pruneInfo = pruneInfo), spec))
        
        result[['fitML']] <- fitML2
        
      } else {
        message <- "The MCMC-fit found a better likelihood outside of the search-region of the ML-fit."
        cat(message)
        warning(message)
      }
    }
    
    result[['fitMCMC']] <- fitMCMC
  }
  class(result) <- c('POUMM', class(result))
  result
}


#' Extract maximum likelihood and degrees of freedom from a fitted POUMM model
#' @param object An object of class POUMM.
#' @param ... not used; included for compliance with generic function logLik.
#' @export
logLik.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    lik <- object$fitML$value
    attr(lik, "df") <- length(object$spec$parLower)
    attr(lik, "nobs") <- object$N
    class(lik) <- "logLik"
    lik
  } else {
    stop("logLik.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood fitted parameters (coefficients) from a fitted 
#' POUMM  model.
#' 
#' @param object An object of class POUMM.
#' @param mapped Logical indicating whether the standard POUMM parameters should
#' also be extracted.
#' @param ... Not used; added for compatibility with generic function coef.
#' @details The parameters extracted are the ones used as input to the model's
#' parMapping function. 
#' 
#' 
#' @return A named vector with the fitted parameters of the model.
#' 
#' @importFrom stats coef
#' 
#' @export
coef.POUMM <- function(object, mapped = FALSE, ...) {
  if("POUMM" %in% class(object)) {
    pars <- object$fitML$par
    g0 <- attr(object$fitML$value, "g0")
    if(mapped) {
      pars <- object$spec$parMapping(pars)
    }
    if("g0" %in% names(pars)) {
      pars["g0"] <- g0
    }
      
    pars
  } else {
    stop("coef.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood expected genotypic values at the tips of a tree,
#' to which a POUMM model has been previously fitted
#' @param object An object of class POUMM.
#' @param vCov A logical indicating whether a list with the genotypic values and their variance covariance matrix should be returned or only a vector of the genotypic values (default is FALSE).
#' @param ... Not used; added for compatibility with generic function fitted.
#' @return If vCov == TRUE, a list with elements g - the genotypic values and
#' vCov - the variance-covariance matrix of these values for the specific tree,
#' observed values z and POUMM ML-fit. If vCov == FALSE, only the vector of 
#' genotypic values corresponding to the tip-labels in the tree is returned.
#' 
#' @importFrom stats fitted
#' @export
fitted.POUMM <- function(object, vCov=FALSE, ...) {
  if("POUMM" %in% class(object)) {
    g0 <- coef.POUMM(object, mapped=TRUE)['g0']
    if(is.nan(g0)) {
      warning("Genotypic values cannot be inferred for g0=NaN; Read documentaton for parMapping in ?specifyPOUMM and use a parMapping function that sets a finite value or NA for g0.")
    }
    p <- coef(object, mapped = TRUE)
    gList <- gPOUMM(object$pruneInfo$z, object$pruneInfo$tree, g0,
                    p["alpha"], p["theta"], p["sigma"], p["sigmae"])
    g = as.vector(gList$mu.g.poumm)
    names(g) <- object$pruneInfo$tree$tip.label
    if(vCov) {
      list(g = g, vCov = gList$V.g.poumm)
    } else {
      g
    }
  } else {
    stop("fitted.POUMM called on non POUMM-object.")
  }
}

#' Number of tips in a phylogenetic tree, POUMM has been fit on.
#' @param object An object of class POUMM.
#' @param ... Not used; added for compatibility with generic function nobs.
#' @return The number of tips in the tree, POUMM has been called on
#' 
#' @export
nobs.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    object$N
  } else {
    stop("nobs.POUMM called on a non POUMM-object.")
  }  
}

#' Extract maximum likelihood environmental contributions (residuals) at the tips of a tree, to which a POUMM model has been fitted.
#' @param object An object of class POUMM.
#' @param ... Not used; added for compatibility with generic function residuals.
#' @return The vector of e-values (residuals) corresponding to the tip-labels in
#'  the tree.
#'  
#' @importFrom stats fitted
#' 
#' @export
residuals.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    e <- object$pruneInfo$z - fitted(object)
    names(e) <- object$pruneInfo$tree$tip.label
    e
  } else {
    stop("residuals.POUMM called on a non POUMM-object.")
  }
}

#' Plots of a POUMM-fit
#' @param x An object of class POUMM.
#' @param type A character indicating the type of plot(s) to be generated.
#'   Defaults to "MCMC", resulting in a trace and density plot for the selected
#'   statistics (see argument stat).
#' @param doPlot Logical indicating whether a plot should be printed on the 
#'   currently active graphics device or whether to return a list of ggplot 
#'   objects for further processing. Defaults to TRUE.
#' @param interactive Logical indicating whether the user should press a key 
#'   before generating a next plot (when needed to display two or more plots).
#'   Defaults to TRUE. Meaningless if 
#'   doPlot = FALSE.
#' @param stat A character vector with the names of statistics to be plotted.
#'   These should be names from the stats-list (see argument statFunctions).
#'   Defaults to c("alpha", "theta", "sigma", "sigmae", "H2tMean", "H2tInf").
#' @param chain A vector of integers indicating the chains to be plotted. 
#' @param startMCMC,endMCMC,thinMCMC Integers used to extract a sample from the 
#'   MCMC-chain; passed to summary().
#' @param statFunctions Named list of statistics functions; passed to summary().
#' @param doZoomIn (type MCMC only) A logical value indicating whether the 
#'   produced plots should have a limitation on the x-axis according to an 
#'   expression set in zoomInFilter (see below). Default value is FALSE.
#' @param zoomInFilter A character string which evaluates as logical value. If 
#'   doZoomIn is set to TRUE, this filter is applied to each point in each MCMC
#'   chain and the data-point is filtered out if it evaluates to FALSE. This 
#'   allows to zoomIn the x-axis of density plots but should be use with caution,
#'   since filtering out points from the MCMC-sample can affect the kernel 
#'   densities. Default value is 
#'    paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
#'           " (value >= HPDLower & value <= HPDUpper))"). 
#'   The identifiers in this expression can be any
#'   column names found in a summary of a POUMM object.
#' @param ... not used, needed for consistency with the generic plot-function.
#'
#' @return If doPlot==FALSE, a named list containing a member called data of
#'   class data.table and several members of class ggplot.
#' @export
plot.POUMM <- 
  function(x, type=c("MCMC"), 
           doPlot = TRUE, interactive = TRUE,
           stat=c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
           chain=NULL,
           startMCMC = NA, endMCMC = NA, thinMCMC = 1000, 
           statFunctions = statistics(x),
           doZoomIn = FALSE,
           zoomInFilter = 
             paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
                    " (value >= HPDLower & value <= HPDUpper))"),
           ...) {
  
  if("POUMM" %in% class(x)) {
    summ <- summary(x, mode = "expert", 
                    startMCMC = startMCMC, endMCMC = endMCMC, 
                    thinMCMC = thinMCMC, stats = statFunctions)
    plot(summ, doPlot = doPlot, interactive = interactive, 
         stat = stat, chain = chain, doZoomIn = doZoomIn, zoomInFilter = zoomInFilter)
  } else {
    stop("plot.POUMM called on a non POUMM-object.")
  }
}