\name{affycompPlot}
\alias{affycompPlot}
\alias{affycomp.figure.calls}
\alias{affycomp.compfigs.calls}
\alias{affycomp.figures}
\alias{affycomp.compfigs}
\title{Assessment Plots}
\description{
Function that makes assessment plot}
}
\usage{
affycompPlot(...,assessment.list=NULL,method.names=NULL,
             figure1.xlim=c(-4,15),figure1.ylim=c(-10,12),
             figure1b.xlim=c(-2,14),figure1b.ylim=c(-6,5),
             figure6a.xlim=c(-12,12),figure6a.ylim=c(-12,12),
             figure6b.xlim=c(-3,3),figure6b.ylim=c(-6,6))

affycomp.compfigs(l, method.names = NULL, figure1.xlim = c(-4, 15),
                 figure1.ylim = c(-10, 12), figure1b.xlim = c(-4, 15),
                 figure1b.ylim = c(-4, 4), figure6a.xlim = c(-12, 12),
                 figure6a.ylim = c(-12, 12), figure6b.xlim = c(-3, 3),
                 figure6b.ylim = c(-6, 6))
affycomp.figures(l)
affycomp.figure.calls(what)
affycomp.compfigs.calls(what)
}
\arguments{
    \item{\dots}{lists produced by the assessment functions (one for
    each method) separated by commas.}
  \item{assessment.list}{Alternatively, one can also send a list of
    lists produced by one of the assessment functions}.
  \item{method.names}{A character vector with the names of the
    epxression measure methodology.}
  \item{figure1.xlim}{x-axis lim used for the plots in Figure 1.}
  \item{figure1.ylim}{y-axis lim used for the plots in Figure 1.}
  \item{figure1b.xlim}{x-axis lim used for the plots in Figure 1b.}
  \item{figure1b.ylim}{y-axis lim used for the plots in Figure 1b.}
  \item{figure6a.xlim}{x-axis lim used for the plots in Figure 6a.}
  \item{figure6a.ylim}{y-axis lim used for the plots in Figure 6a.}
  \item{figure6b.xlim}{x-axis lim used for the plots in Figure 6b.}
  \item{figure6b.ylim}{y-axis lim used for the plots in Figure 6b.}
  \item{l}{list with assessment lists as components.}
  \item{what}{a dummy variable to know what function call to create.}
}

\details{
Read the vignette for more details on what each Figure is. Once an
assessment is used this function knows what to do. You can call any of
the assessment functions described in \code{\link{assessSpikeIn}},
\code{\link{assessDilution}} and \code{\link{assessSD}}.

\code{affycomp.figures}, \code{affycomp.figure.calls},
\code{affycomp.compfigs.calls} are auxiliary functions.
}
\value{
  Figures are produced.
}
\author{Rafael A. Irizarry}
\examples{
data(rma.assessment)
data(mas5.assessment)
affycompPlot(rma.assessment,mas5.assessment)
affycompPlot(rma.assessment$FC)
affycompPlot(rma.assessment$Dilution,mas5.assessment$Dilution)
affycompPlot(rma.assessment$Dilution,mas5.assessment$Dilution)
affycompPlot(rma.assessment$Signal,mas5.assessment$Signal)
affycompPlot(rma.assessment$Dilution,mas5.assessment$Dilution)
affycompPlot(rma.assessment$FC2,mas5.assessment$FC2)
}
\keyword{hplot}


