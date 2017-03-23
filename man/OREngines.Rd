\name{OREngines}
\docType{data}
\alias{OREngines}
\title{Maintenance times and types of off road engines}
\description{
	This data set refers to 141 diesel engines of trucks used in Brazilian mining. Failure times and planned PM times are recorded in working hours. 208 failure times and 52 PM times were recorded. The data set does not exactly corresponds to the one used by M.L.G. de Toledo et al.. In fact, in their paper the authors have considered AGAN PM, then they have transformed the data set corresponding to the 143 engines with 50 PM into an equivalent data base of 193=141+52 systems, using the fact that after each PM the engine is equivalent to a new one. But when making this transformation they have made some little typos. But, this data set corresponds to the original and consequently correct data set. 
}
\usage{OREngines}
\format{A data frame with 260 lines and 3 columns. Each line corresponds to an event. For each event, the columns \code{System}, \code{Time} and \code{Type} represent respectively the corresponding system (an integer between 1 and 141), the times and the type (-1 for a CM, 1 for a PM).}
\references{
  M.L.G. de Toledo, M.A. Freitas, E.A. Colosimo  and L.G. Gilardoni (2016),
  "Optimal periodic maintenance policy under imperfect repair: A case study on the engines of off-road vehicles",\emph{IIE Transactions} 48(8), pp. 747-758.
}
\seealso{
\code{\link{OREngines_CM}} for the data set in which PM are considered as renewing.
\code{\link{OREngines2_CM}} for the data set studied by M.L.G. de Toledo et al..
}
\keyword{datasets}