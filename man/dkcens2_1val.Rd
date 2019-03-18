\name{dkcens2_1val}
\alias{dkcens2_1val}
\title{
Density evaluation once
}
\description{
Computes the likelihood contribution for one data point in the case of censoring.
}
\usage{
dkcens2_1val(xleft, xright, c_code, distr = NULL, mu = NULL, sigma = NULL)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (xleft, xright, c_code, distr = NULL, mu = NULL, sigma = NULL) 
{
    if (c_code == 1) 
        dk(x = xleft, distr, mu, sigma)
    else if (c_code == 2) 
        pk(xright, distr, mu, sigma)
    else if (c_code == 0) 
        1 - pk(xleft, distr, mu, sigma)
    else if (c_code == 3) 
        pk(xright, distr, mu, sigma) - pk(xleft, distr, mu, sigma)
    else stop("Wrong integer code for censored data")
  }
}
\keyword{internal}
