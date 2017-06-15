## ----setup, include = FALSE----------------------------------------------
# Make results reproducible
set.seed(1)
library(ggplot2)
library(data.table)
knitr::opts_chunk$set(cache = FALSE)
options(digits = 4)

# set this to FALSE to disable cache and run all MCMC-fits.
useCachedResults <- file.exists("UserGuideCache.RData") & TRUE


## ----install, eval=FALSE-------------------------------------------------
#  install.packages('POUMM')
#  install.packages("TreeSim")
#  install.packages("data.table")
#  install.packages("ggplot2")
#  install.packages("lmtest")

## ------------------------------------------------------------------------
N <- 500
g0 <- 0           
alpha <- .5        
theta <- 2        
sigma <- 0.2     
sigmae <- 0.2 

## ---- include=FALSE, eval=useCachedResults-------------------------------
load("UserGuideCache.RData")

## ---- echo=FALSE, fig.height=4.6, fig.width=7, fig.cap="Dashed black and magenta lines denote the deterministic trend towards the long-term mean $\\theta$, fixing the stochastic parameter $\\sigma=0$."----
tStep <- 0.025
t <- seq(0, 6, by = tStep)

plot(t, POUMM::rTrajectoryOU(g0, tStep, alpha, theta, sigma, length(t)), type = 'l', main = "Random OU trajectories", ylab = "g", ylim = c(0, 4))
lines(t, POUMM::rTrajectoryOU(g0, tStep, alpha, theta, 0, length(t)), lty = 2)

lines(t, POUMM::rTrajectoryOU(g0, tStep, alpha*2, theta, sigma, length(t)), col = "magenta")
lines(t, POUMM::rTrajectoryOU(g0, tStep, alpha*2, theta, 0, length(t)), lty = 2, col = "magenta")

lines(t, POUMM::rTrajectoryOU(g0, tStep, alpha, theta, sigma*2, length(t)), col = "blue")

abline(h=theta, lty = 3, col = "darkgrey")

legend("topleft", 
       legend = c(expression(list(alpha == .5, theta == 2, sigma == 0.2)),
                  expression(list(alpha == .5, theta == 2, sigma == 0.4)),
                  expression(list(alpha == .5, theta == 2, sigma == 0)),
                  
                  expression(list(alpha == 1, theta == 2, sigma == 0.2)),
                  expression(list(alpha == 1, theta == 2, sigma == 0)),
                  
                  expression(theta == 2)),
       lty = c(1, 1, 2, 1, 2, 3), 
       col = c("black", "blue", "black", "magenta", "magenta", "darkgrey"))

## ----simulate-tree, results="hide", eval=!useCachedResults---------------
#  # Number of tips
#  tree <- TreeSim::sim.bdsky.stt(N, lambdasky = 1.6, deathsky = .6,
#                                 timesky=c(0, Inf), sampprobsky = 1)[[1]]

## ----simulate-gez-OU, eval=!useCachedResults-----------------------------
#  # genotypic (heritable) values
#  g <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma)
#  
#  # environmental contributions
#  e <- rnorm(length(g), 0, sigmae)
#  
#  # phenotypic values
#  z <- g + e

## ----violin-plots, fig.show = "hold", fig.height=4, fig.width=7, fig.cap="Distributions of the trait-values grouped according to their root-tip distances."----
# This is easily done using the nodeTimes utility function in combination with
# the cut-function from the base package.
data <- data.table(z = z[1:N], t = POUMM::nodeTimes(tree, tipsOnly = TRUE))
data <- data[, group := cut(t, breaks = 5, include.lowest = TRUE)]

ggplot(data = data, aes(x = t, y = z, group = group)) + 
  geom_violin(aes(col = group)) + geom_point(aes(col = group), size=.5)

## ----MaintainCache, echo=FALSE, warning=FALSE, results="hide", message=FALSE, eval=TRUE----
if(!useCachedResults) {
  # Perform the heavy fits locally. 
  # set up a parallel cluster on the local computer for parallel MCMC:
  cluster <- parallel::makeCluster(parallel::detectCores(logical = FALSE))
  doParallel::registerDoParallel(cluster)
  
  fitPOUMM <- POUMM::POUMM(z[1:N], tree, spec=list(thinMCMC = 1000, parallelMCMC=TRUE), verbose = TRUE)
  fitPOUMM2 <- POUMM::POUMM(z[1:N], tree, spec=list(nSamplesMCMC = 4e5, thinMCMC = 1000, parallelMCMC=TRUE))
  
  specH2tMean <- POUMM::specifyPOUMM_ATH2tMeanSeG0(z[1:N], tree, 
                                                   nSamplesMCMC = 4e5, 
                                                   thinMCMC = 1000, parallelMCMC=TRUE)
  fitH2tMean <- POUMM::POUMM(z[1:N], tree, spec = specH2tMean)
  
  
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.
  parallel::stopCluster(cluster)
  
  # delete some of the slots to make cach-file smaller :-(.
  fitPOUMM$pruneInfo <- fitPOUMM2$pruneInfo <- fitH2tMean$pruneInfo <- NULL
  fitPOUMM$loglik <- fitPOUMM2$loglik <- fitH2tMean$loglik <- NULL
  fitPOUMM$fitMCMC$post <- fitPOUMM2$fitMCMC$post <- fitH2tMean$fitMCMC$post <- NULL
  fitPOUMM$spec$parPriorMCMC <- fitPOUMM2$spec$parPriorMCMC <- fitH2tMean$spec$parPriorMCMC <- NULL
  fitPOUMM$spec$parInitMCMC <- fitPOUMM2$spec$parInitMCMC <- fitH2tMean$spec$parInitMCMC <- NULL
  
  save(g, z, tree, e, 
       fitPOUMM2, fitPOUMM, fitH2tMean,
       file="UserGuideCache.RData")
} 
# restore the pruneInfo since it is needed afterwards.
fitPOUMM$pruneInfo <- fitPOUMM2$pruneInfo <- fitH2tMean$pruneInfo <- POUMM::pruneTree(tree, z[1:N])

## ----fitPOUMM-1, results="hide", message=FALSE, warning=FALSE, eval=FALSE----
#  fitPOUMM <- POUMM::POUMM(z[1:N], tree)

## ---- fig.height=5.4, fig.show="hold", fig.width=7.2, warning=FALSE, fig.cap="MCMC traces and posterior density plots from a POUMM MCMC-fit. Black dots on the x-axis indicate the ML-fit.", results="hide", eval=TRUE----
plot(fitPOUMM, showUnivarDensityOnDiag = TRUE)

## ---- warning=FALSE, eval=TRUE-------------------------------------------
summary(fitPOUMM)

## ----fitPOUMM-2, results="hide", eval=FALSE------------------------------
#  fitPOUMM2 <- POUMM::POUMM(z[1:N], tree, spec=list(nSamplesMCMC = 4e5))

## ---- fig.height=5.4, fig.show="hold", fig.width=7.2, warning=FALSE, results="hide", eval=TRUE----
pl <- plot(fitPOUMM2, doPlot = FALSE, doZoomIn = TRUE)
pl$densplot

## ---- warning=FALSE, eval=TRUE-------------------------------------------
summary(fitPOUMM2)

## ------------------------------------------------------------------------
tMean <- mean(POUMM::nodeTimes(tree, tipsOnly = TRUE))
tMax <- max(POUMM::nodeTimes(tree, tipsOnly = TRUE))

c(# phylogenetic heritability at mean root-tip distance: 
  H2tMean = POUMM::H2(alpha, sigma, sigmae, t = tMean),
  # phylogenetic heritability at long term equilibirium:
  H2tInf = POUMM::H2(alpha, sigma, sigmae, t = Inf),
  # empirical (time-independent) phylogenetic heritability, 
  H2e = POUMM::H2e(z[1:N], sigmae),
  # genotypic variance at mean root-tip distance: 
  sigmaG2tMean = POUMM::varOU(t = tMean, alpha, sigma),
  # genotypic variance at max root-tip distance: 
  sigmaG2tMean = POUMM::varOU(t = tMax, alpha, sigma),
  # genotypic variance at long-term equilibrium:
  sigmaG2tInf = POUMM::varOU(t = Inf, alpha, sigma)
  )

## ---- warning=FALSE------------------------------------------------------
c(H2empirical = var(g[1:N])/var(z[1:N]))
summary(fitPOUMM2)["H2tMean"==stat, unlist(HPD)]

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # set up a parallel cluster on the local computer for parallel MCMC:
#  cluster <- parallel::makeCluster(parallel::detectCores(logical = FALSE))
#  doParallel::registerDoParallel(cluster)
#  
#  fitPOUMM <- POUMM::POUMM(z[1:N], tree, spec=list(parallelMCMC = TRUE))
#  
#  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.
#  parallel::stopCluster(cluster)

## ------------------------------------------------------------------------
specPMM <- POUMM::specifyPMM(z[1:N], tree)
fitPMM <- POUMM::POUMM(z[1:N], tree, spec = specPMM, doMCMC=FALSE)

## ------------------------------------------------------------------------
lmtest::lrtest(fitPMM, fitPOUMM2)

## ------------------------------------------------------------------------
gBM <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha = 0, theta = 0, sigma = sigma)
zBM <- gBM + e

fitPMM_on_zBM <- POUMM::POUMM(zBM[1:N], tree, spec = specPMM, doMCMC = FALSE)
fitPOUMM_on_zBM <- POUMM::POUMM(zBM[1:N], tree, doMCMC = FALSE)

lmtest::lrtest(fitPMM_on_zBM, fitPOUMM_on_zBM)

## ---- message=FALSE, warning=FALSE, eval=TRUE----------------------------
specH2tMean <- POUMM::specifyPOUMM_ATH2tMeanSeG0(z[1:N], tree, nSamplesMCMC = 4e5)
# Mapping from the sampled parameters to the standard POUMM parameters:
specH2tMean$parMapping
# Prior for the MCMC sampling
specH2tMean$parPriorMCMC
# Bounds for the maximum likelihood search
specH2tMean$parLower
specH2tMean$parUpper

## ----eval=FALSE----------------------------------------------------------
#  fitH2tMean <- POUMM::POUMM(z[1:N], tree, spec = specH2tMean)

## ---- fig.height=5.4, fig.show="hold", fig.width=7.2, warning=FALSE------
plot(fitH2tMean, stat = c("H2tMean", "H2e", "H2tInf", "sigmae"), 
     showUnivarDensityOnDiag = TRUE, 
     doZoomIn = TRUE, doPlot = TRUE)

## ---- warning=FALSE------------------------------------------------------
summary(fitH2tMean)[stat %in% c("H2tMean", "H2e", "H2tInf", "sigmae")]

## ------------------------------------------------------------------------
# Compare global empirical heritability
H2eGlobal <- POUMM::H2e(z[1:N], sigmae = coef(fitH2tMean)['sigmae'])
# versus recent empirical heritability
H2eRecent <- POUMM::H2e(z[1:N], tree, sigmae = coef(fitH2tMean)['sigmae'], tFrom = 5)
print(c(H2eGlobal, H2eRecent))

## ----create-references, echo=FALSE, include=FALSE, eval=TRUE-------------
likCalculation <- c("Rcpp", "RcppArmadillo", "Rmpfr")
mcmcSampling <- c("adaptMCMC")
mcmcDiagnosis <- c("coda")
otherPackages <- c("parallel", "foreach", "data.table", "Matrix", "gsl")
treeProcessing <- c("ape")
reporting <- c("data.table", "ggplot2", "GGally", "lmtest")
testing <- c("testthat", "mvtnorm", "TreeSim")
 
packagesUsed <- c(likCalculation, mcmcDiagnosis, otherPackages, treeProcessing, reporting, testing)

printPackages <- function(packs) {
  res <- ""
  for(i in 1:length(packs)) {
    res <- paste0(res, paste0(packs[i], ' v', packageVersion(packs[i]), ' [@R-', packs[i], ']'))
    if(i < length(packs)) {
      res <- paste0(res, ', ')
    }
  }
  res
}

# Write bib information (this line is executed manually and the bib-file is edited manually after that)
# knitr::write_bib(packagesUsed, file = "./REFERENCES-R.bib")

