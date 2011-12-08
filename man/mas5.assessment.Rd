\name{mas5.assessment}
\docType{data}
\alias{mas5.assessment}

\title{Example of the result of assessments}

\description{
  The Dilution and both (HGU95 and HGU133) types of Spike-in
  data were processed with Affymetrix MAS 5.0 software,
  yielding three "MAS 5.0" \link[Biobase]{ExpressionSet}'s.
  (These are available, in csv-format, at
  \url{http://affycomp.jhsph.edu/AFFY2/rafa@jhu.edu/030424.1033/}.)
  Then various assessment functions from the affycomp package
  (most recently, version 1.28.0) were applied.
  \code{mas5.assessment} resulted from
  \code{\link[affycomp]{assessAll}} on Dilution and HGU95.
  See \link[affycompData]{mas5.assessment} in affycompData
  for results of other assessments.
}

\usage{data(mas5.assessment)
}
\format{A list of list.}


\keyword{datasets}
