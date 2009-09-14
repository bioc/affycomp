\name{assessSpikeIn}
\alias{assessSpikeIn}
\alias{assessFC}
\alias{assessFC2}
\alias{assessMA}
\alias{assessSignal}

\title{Assessment functions for Spike In Data}
\description{
  These functions are assessment functions.
  Each takes an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} and returns a list with necessary
  information to create assessment plots.}
\usage{
assessSpikeIn(s,method.name=NULL,verbose=TRUE)

assessMA(exprset,method.name=NULL)
assessSignal(exprset,method.name=NULL)
assessFC(exprset,method.name=NULL)
assessFC2(exprset,method.name=NULL)
}
\arguments{
    \item{s}{An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
      containing 
     expression measures for Affymetrix's spike-in experiment.}
  \item{exprset}{An
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing 
     expression measures for Affymetrix's spike-in experiment.}
  \item{method.name}{Name of expression measure being assessed.}
  \item{verbose}{logical. If \code{TRUE} show messages.}
}
\details{
  \code{assessMA} performs the assessment for the MA-plot (Figure 1),
  \code{assessSignal} performs the assessment for signal detection plot
  (Figure 4a), \code{assessFC} performs assessments used by fold-change
  related  plots (Figures 5a, 6a, 6b). \code{assessFC2} is for the ROC
for genes with nominal fold changes of 2 (Figure
5b). \code{assessSpikeIn} is a wrapper for all these and returns a list
of lists. }
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}
\keyword{manip}
