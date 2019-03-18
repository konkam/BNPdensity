\name{qt_}
\alias{qt_}
\title{
Quantile function non-standard Student-t}
\description{
Compuytes the quantiles.
}
\usage{
qt_(p, df, mean, sd)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (p, df, mean, sd) 
{
    sd * qt(p, df, ncp = 0) + mean
  }
}
\keyword{internal}
