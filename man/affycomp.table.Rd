\name{affycompTable}
\alias{affycompTable}
\alias{tableAll}
\alias{tableDilution}
\alias{tableFC}
\alias{tableFC2}
\alias{tableSignal}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Expression Assessment Table}
\description{
  These functions takes as an argument the output of the asessment
  functions.
}
\usage{
affycompTable(...,Table=NULL,assessment.list=NULL,method.names=NULL)

tableAll(...,assessment.list=NULL,method.names=NULL)

tableDilution(l, method.names=NULL)

tableFC(l, method.names=NULL)

tableFC2(l, method.names=NULL)

tableSignal(l, method.names=NULL)
}
\arguments{
  \item{\dots}{lists produced by the assessment functions}
  \item{Table}{If \code{TableAll} was used one can send it through this argument}
  \item{assessment.list}{Alternatively, one can also send a list of lists produced by \code{\link{tableAll}}.}
  \item{method.names}{A character vector with the names of the
    epxression measure methodology.}
  \item{l}{list of assessments.}
}
\details{
  Read the vignette for more details on what the entries of the table
  are. \code{affycompTable} has a few entries per graph. \code{tableAll}
  has more entries. Once an
  assessment is used this function knows what to do. You can call any of
  the assessment functions described in \code{\link{assessSpikeIn}},
  \code{\link{assessDilution}} and \code{\link{assessSD}}.
  
}
\value{
A matrix. One column per each method and one row for each comparison.}
}
\author{Rafael A. Irizarry}
\examples{
data(rma.assessment) ##this was produced with affycomp.assess
data(mas5.assessment) ##this one too
tmp <- affycompTable(mas5.assessment,rma.assessment)
format(tmp,digit=2)
}
\keyword{manip}


