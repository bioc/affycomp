\name{assessAll}
\alias{assessAll}
\alias{affycomp}
\title{Assessment functions}
\description{
  Assessment functions.
  Takes a couple of 
  \code{\link[Biobase]{exprSet-class}}, one for spike in another for the
  dilution and returns a list with necessary
  information to create assessment plots.}
\usage{
assessAll(d,s,method.name=NULL,verbose=TRUE)

affycomp(d,s,method.name=NULL,verbose=TRUE,return.it=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
    \item{d}{An \code{\link[Biobase]{exprSet-class}} containing
    the exrpession measures for the Gene Logic's dilution experiment.}
  \item{s}{An \code{\link[Biobase]{exprSet-class}} containing
    expression measures for the 
    Affymetrix's spike-in experiment.}
  \item{method.name}{Name of expression measure being assessed.}
  \item{verbose}{verbosity flag.}
\item{return.it}{if \code{TRUE} returns assessment lists.}  
}
\details{
  \code{assessAll} performs assessments for Figures 1-6. It is a wrapper
  for \code{assessDilution} and \code{assessSpikeIn}.

\code{affycomp} is a wrapper that does it all... including the plottitng
and creation of table.
}
\value{
  Lists with the necessary information to make the Figures.
}
\author{Rafael A. Irizarry}

\keyword{manip}


