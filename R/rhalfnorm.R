rhalfnorm <-
function (n, mean = 0, sd = 1) 
{
    abs(rnorm(n, mean, sd))
}
