\name{readin}
\alias{readin}
\alias{read.spikein}
\alias{read.dilution}
\title{Read Expression Date Sets}
\description{
Reads a comma delimeted file containing the expression values of the
dilution and spike-in data sets and creates a \code{\link[Biobase]{exprSet}}
}
\usage{
read.dilution(filename)
read.spikein(filename)
}
\arguments{
  \item{filename}{character containing the filename to be read.}
}
\details{
 The file to be read must be come delimeted with the first row
 containing the cel fileneames (case sensitive). The first column must
 be the Affymetrix gene identifiers. \code{read.dilution} will put
 things in the right place.
}
\value{An \code{\link[Biobase]{exprSet}}.
}
\author{Rafael A. Irizarry}
\keyword{manip}
