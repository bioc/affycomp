\name{readin}
\alias{readin}
\alias{read.spikein}
\alias{read.newspikein}
\alias{read.dilution}
\title{Read Expression Date Sets}
\description{
Reads a comma delimeted file containing the expression values of the
dilution and spike-in data sets and creates a \code{\link[Biobase]{exprSet}}
}
\usage{
read.dilution(filename)
read.spikein(filename,cdfName=c("hgu95a","hgu133a"),remove.xhyb=TRUE)
read.newspikein(filename)
}
\arguments{
  \item{filename}{character containing the filename to be read.}
\item{cdfName}{are we reading data from the hgu95a or hgu133a spike-in experiment?}
\item{remove.xhyb}{logical. If \code{TRUE} possible cross hybridizers
are removed from the HGU133A spikein. See \code{\link{remove.hgu133a.xhyb}}.}
}
\details{
 The file to be read must be come delimeted with the first row
 containing the cel fileneames (case sensitive). The first column must
 be the Affymetrix gene identifiers. \code{read.dilution} will put
 things in the right place.

 \code{read.newspikein} is a wrapper to read results from the hgu133a
spikein experiment.
}
\value{An \code{\link[Biobase]{exprSet}}.
}
\author{Rafael A. Irizarry}
\keyword{manip}
