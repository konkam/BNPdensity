qhalfnorm <-
  function(p, mean = 0, sd = 1) {
    qnorm(
      p * (1 - pnorm(0, mean, sd)) + pnorm(0, mean, sd),
      mean, sd
    )
  }
