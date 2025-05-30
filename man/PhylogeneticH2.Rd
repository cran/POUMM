% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paramsPOUMM.R
\name{PhylogeneticH2}
\alias{PhylogeneticH2}
\alias{alpha}
\alias{sigmaOU}
\alias{sigmae}
\alias{H2e}
\title{Phylogenetic Heritability}
\usage{
alpha(H2, sigma, sigmae, t = Inf)

sigmaOU(H2, alpha, sigmae, t = Inf)

sigmae(H2, alpha, sigma, t = Inf)

H2e(z, sigmae, tree = NULL, tFrom = 0, tTo = Inf)
}
\arguments{
\item{H2}{Phylogenetic heritability at time t.}

\item{sigmae}{Numeric, environmental phenotypic deviation at the tips.}

\item{t}{Numeric value denoting evolutionary time (i.e. distance from the 
root of a phylogenetic tree).}

\item{alpha, sigma}{Numeric values or n-vectors, parameters of the OU process;
alpha and sigma must be non-negative. A zero alpha is interpreted as the 
Brownian motion process in the limit alpha -> 0.}

\item{z}{Numerical vector of observed phenotypes.}

\item{tree}{A phylo object.}

\item{tFrom, tTo}{Numerical minimal and maximal root-tip distance to limit the
calculation.}
}
\value{
All functions return numerical values or NA, in case of invalid 
  parameters
}
\description{
The phylogenetic heritability, \eqn{H^2}, is defined as the 
  ratio of the genetic variance over the total phenotypic variance expected 
  at a given evolutionary time t (measured from the root of the tree). Thus,
  the phylogenetic heritability connects the parameters alpha, sigma and
  sigmae of the POUMM model through a set of equations. The functions
  described here provide an R-implementation of these equations.
}
\details{
The function sigmae uses the formula H2 = varOU(t, alpha, sigma) /
  (varOU(t, alpha, sigma) + sigmae^2)
}
\section{Functions}{
\itemize{
\item \code{alpha()}: Calculate alpha given time t, H2, sigma and sigmae

\item \code{sigmaOU()}: Calculate sigma given time t, H2 at time t, alpha
and sigmae

\item \code{sigmae()}: Calculate sigmae given alpha, sigma, and H2 at
time t

\item \code{H2e()}: "Empirical" phylogenetic heritability estimated
from the empirical variance of the observed phenotypes and sigmae

}}
\note{
This function is called sigmaOU and not simply sigma to avoid a conflict 
with a function sigma in the base R-package.
}
\examples{
# At POUMM stationary state (equilibrium, t=Inf)
H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = Inf)     # 0.4
alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = Inf)    # 0.75
sigma <- sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = Inf) # 1
sigmae <- sigmae(H2 = H2, alpha = 0.75, sigma = 1, t = Inf) # 1

# At finite time t = 0.2
H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = 0.2)     # 0.1473309
alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = 0.2)    # 0.75
sigma <- sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = 0.2) # 1
sigmae <- sigmae(H2  =  H2, alpha = 0.75, sigma = 1, t = 0.2) # 1

   
}
\seealso{
OU
}
