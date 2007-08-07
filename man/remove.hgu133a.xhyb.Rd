\name{remove.hgu133a.xhyb}
\alias{remove.hgu133a.xhyb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove crosshybridizers}
\description{
This functions removes possible cross hybridizers from Affymetrix
HGU133A spike-in experiment
}
\usage{
remove.hgu133a.xhyb(s, bp = c("200", "150", "100"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{an \code{exprSet} containing the HGU133A spike-in}
  \item{bp}{number of base pair matches needed to define a possible
    cross hybridizer}
}
\details{Some details are contained in the help file for \code{\link{hgu133a.spikein.xhyb}}}

\value{
  An \code{exprSet} with probeset removed}
\seealso{ \code{\link{hgu133a.spikein.xhyb}}}
\keyword{manip}
