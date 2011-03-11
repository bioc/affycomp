\name{assessDilution}
\alias{assessDilution}

\title{Assessment functions for Dilution Data}
\description{
  Assessment function. Takes an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} and returns a
  list with necessary information to create assessment plots.}
\usage{
assessDilution(exprset,method.name=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exprset}{An
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing 
     expression measures for GeneLogic's dilution experiment.}
  \item{method.name}{Name of expression measure being assessed.}
}
\details{
  \code{assessDilution} performs the assessment for the plots related to
  Dilution (Figures 2, 3, 4b)
}
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}
\keyword{manip}
