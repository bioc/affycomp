\name{assessDilution}
\alias{assessDilution}


\title{Assessment functions for Dilution Data}
\description{
  Assessment functoion.
  Takes an
  \code{\link[Biobase]{exprSet-class}} and returns a list with necessary
  information to create assessment plots.}
\usage{
assessDilution(exprset,method.name=NULL)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exprset}{An \code{\link[Biobase]{exprSet-class}} containing
     exrpession measures for 
    Affymetrix's spike-in experiment.}
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


