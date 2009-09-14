\name{assessSD}
\alias{assessSD}

\title{SD Assessment functions}
\description{
  Assessment function for standard deviation estimates.
  Takes a dilution data
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} and returns a
  list with necessary information to create assessment plot.}
\usage{
assessSD(exprset, method.name=NULL, logx=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exprset}{An
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing 
     expression measures for Affymetrix's spike-in experiment.}
  \item{method.name}{Name of expression measure being assessed.}
  \item{logx}{Logical. If \code{TRUE} expression is logged for plot. See
    details.} 
}
\details{
  \code{assessSD} does the assessment for Figure 7. This requires the
  \code{ExpressionSet} to have standard error estimates for the expression
  measure. Some expression (such as dChip) will have SEs for the
  original scale. Others, like RMA will have them for the log scale. For
  original scales, making \code{logx=TRUE} is recommended.
}
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}
\keyword{manip}
