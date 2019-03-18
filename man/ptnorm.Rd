\name{ptnorm}
\alias{ptnorm}
\title{
Distribution function truncated normal
}
\description{
Computes the cumulative distribution function.
}
\usage{
ptnorm(q, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, log.p = FALSE)
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
function (q, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, 
    log.p = FALSE) 
{
    ret <- numeric(length(q))
    if (lower.tail) {
        ret[q < lower] <- 0
        ret[q > upper] <- 1
    }
    else {
        ret[q < lower] <- 1
        ret[q > upper] <- 0
    }
    ret[upper < lower] <- NaN
    ind <- q >= lower & q <= upper
    if (any(ind)) {
        denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, 
            sd)
        if (lower.tail) 
            qtmp <- pnorm(q, mean, sd) - pnorm(lower, mean, sd)
        else qtmp <- pnorm(upper, mean, sd) - pnorm(q, mean, 
            sd)
        if (log.p) 
            qtmp <- log(qtmp) - log(denom)
        else qtmp <- qtmp/denom
        ret[q >= lower & q <= upper] <- qtmp[ind]
    }
    ret
  }
}
\keyword{internal}
