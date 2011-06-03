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
     expression measures for Affymetrix's dilution experiment.}
  \item{method.name}{Name of expression measure being assessed.}
  \item{logx}{
  Logical indicating whether the average expression being computed should
  be logged, default no. See details.
  } 
}
\details{
  \code{assessSD} does the assessment for Figure 7. This requires the
  \code{ExpressionSet} to have standard error estimates for the expression
  measure. Some expression measures (e.g. dChip) will have SEs in original
  scale, where others (e.g. RMA) will have them in log scale. For the former,
  use \code{logx=TRUE}.
}
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}
\keyword{manip}
