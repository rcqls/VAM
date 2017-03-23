\name{OREngines_CM}
\docType{data}
\alias{OREngines_CM}
\title{Maintenance times and types of off road engines (with PM renewing assumption)}
\description{
	This data set refers to 141 diesel engines of trucks used in Brazilian mining. Failure times and planned PM times are recorded in working hours. 208 failure times and 52 PM times were recorded. Similarly to what has been done by M.L.G. de Toledo et al., PM are supposed to renew the corresponding system. The data set is then transformed in a sample of 191 sequences of event times (191=141+52-2: 141 systems, 52 AGAN PM and 2 systems with only 1 PM). This data set does not exactly corresponds to the one used by M.L.G. de Toledo et al.. In fact, in their paper the authors have made some little typos when they have transformed the original data set. This data set corresponds to the correct data set transformation. 
}
\usage{OREngines_CM}
\format{A data frame with 191 lines and 3 columns. Each line corresponds to an event. For each event, the columns \code{System}, \code{Time} and \code{Type} represent respectively the corresponding trajectory (an integer between 1 and 191), the times and the type (-1 for a CM).}
\references{
  M.L.G. de Toledo, M.A. Freitas, E.A. Colosimo  and L.G. Gilardoni (2016),
  "Optimal periodic maintenance policy under imperfect repair: A case study on the engines of off-road vehicles",\emph{IIE Transactions} 48(8), pp. 747-758.
}
\seealso{
\code{\link{OREngines}} for the original data set with PM times.
\code{\link{OREngines2_CM}} for the data set studied by M.L.G. de Toledo et al..
}
\keyword{datasets}