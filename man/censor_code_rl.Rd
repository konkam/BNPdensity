\name{censor_code_rl}
\alias{censor_code_rl}
\title{
Censor code right-left
}
\description{
Creates censoring code 0:3.
}
\usage{
censor_code_rl(left, right)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (left, right) 
{
    test_ = function(k) {
        if (is.na(left[[k]]) & is.na(right[[k]])) 
            NA
        else if (is.na(left[[k]])) 
            2
        else if (is.na(right[[k]])) 
            0
        else if (left[[k]] == right[[k]]) 
            1
        else 3
    }
    sapply(seq_along(left), FUN = test_)
  }
}
\keyword{internal}
