\name{qtnorm}
\alias{qtnorm}
\title{
Quantile function truncated normal
}
\description{
Computes the quantiles.
}
\usage{
qtnorm(p, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, log.p = FALSE)
}
\details{
For internal use
}
\references{
Taken from 
}
\author{
C. H. Jackson 
}
\note{
Taken from \code{msm} R-package.
}
\examples{
## The function is currently defined as
function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, 
    log.p = FALSE) 
{
    qgeneric(ptnorm, p = p, mean = mean, sd = sd, lower = lower, 
        upper = upper, lbound = lower, ubound = upper, lower.tail = lower.tail, 
        log.p = log.p)
  }
}
\keyword{internal}
