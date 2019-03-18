qtnorm <-
function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, 
    log.p = FALSE) 
{
    qgeneric(ptnorm, p = p, mean = mean, sd = sd, lower = lower, 
        upper = upper, lbound = lower, ubound = upper, lower.tail = lower.tail, 
        log.p = log.p)
}
