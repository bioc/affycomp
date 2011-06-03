\name{assessSpikeIn2}
\alias{assessSpikeIn2}
\alias{assessSpikeInSD}
\alias{assessLS}
\alias{assessMA2}

\title{New Assessment functions for Spike In Data}
\description{
  These functions are assessment functions.
  Each takes an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} and returns a
  list with necessary 
  information to create assessment plots.}
\usage{
assessSpikeIn2(s, method.name=NULL, verbose=TRUE)
assessSpikeInSD(exprset, method.name=NULL, span=1/3)
assessLS(exprset, method.name=NULL)
assessMA2(exprset, method.name=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
    \item{s}{An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
      containing expression measures for Affymetrix's spike-in experiment.}
  \item{exprset}{An
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing 
     expression measures for Affymetrix's spike-in experiment.}
  \item{method.name}{Name of expression measure being assessed.}
  \item{verbose}{logical. If \code{TRUE} show messages.}
  \item{span}{span used in call to \code{loess}.}
}
\details{
  \code{assessMA2} performs the assessment for the second MA-plot (Figure 1b),
  and \code{assessLS} performs the assessment for signal detection plot
  (Figure 4c).
  \code{assessMA2} also performs assessments used by fold-change related plots
  (Figures 5a,b) and the ROC plots (Figures 5c,d,e).
  \code{assessSpikeInSD} is for the standard deviation assessment in Figure 2b.
  \code{assessSpikeIn2} is a wrapper for all these and returns a list of lists.  
}
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}
\keyword{manip}
