\name{pt_}
\alias{pt_}
\title{
Distribution function non-standard student-t
}
\description{
Computes the cdf.
}
\usage{
pt_(x, df, mean, sd)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (x, df, mean, sd) 
{
    pt((x - mean)/sd, df, ncp = 0)
  }
}
\keyword{internal}
