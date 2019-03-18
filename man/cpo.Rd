\name{cpo}
\alias{cpo}
\title{
Conditional predictive ordinate function}
\description{
This function computes conditional predictive ordinates for each datapoint. 
}
\usage{
cpo(obj)
}
\details{
For internal use. 
}
\examples{
## The function is currently defined as
function (obj) 
{
    fx <- obj$fx
    cpo <- 1/apply(1/fx, 1, mean)
    return(cpo)
  }
}
\keyword{internal}
