\name{OREngines2_CM}
\docType{data}
\alias{OREngines2_CM}
\title{Maintenance times and types of off road engines (with PM renewing assumption)}
\description{
	This data set refers to 141 diesel engines of trucks used in Brazilian mining. Failure times and planned PM times are recorded in working hours. 208 failure times and 52 PM times were recorded. The data set corresponds to the one studied by M.L.G. de Toledo et al. in which PM are supposed to be AGAN. Then, the original data set corresponding to the 143 engines with 50 PM is transformed into an equivalent data base of 193=141+52 systems, using the fact that after each PM the engine is equivalent to a new one. But when making this transformation they have made some little typos. This data set corresponds to the one used by M.L.G. de Toledo et al. and not to the correct data set.
}
\usage{OREngines2_CM}
\format{A data frame with 191 lines and 3 columns. Each line corresponds to an event. For each event, the columns \code{System}, \code{Time} and \code{Type} represent respectively the corresponding trajectory (an integer between 1 and 191), the times and the type (-1 for a CM).}
\references{
  M.L.G. de Toledo, M.A. Freitas, E.A. Colosimo  and L.G. Gilardoni (2016),
  "Optimal periodic maintenance policy under imperfect repair: A case study on the engines of off-road vehicles",\emph{IIE Transactions} 48(8), pp. 747-758.
}
\seealso{
\code{\link{OREngines}} for the original data set with PM times.
\code{\link{OREngines_CM}} for the correct transformed data set.
}
\keyword{datasets}