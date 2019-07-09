phalft <-
  function(q, df = 1, mean = 0, sd = 1) {
    ifelse(q < 0, 0, 1) * (pt_(q, df, mean, sd) - pt_(
      0, df,
      mean, sd
    )) / (1 - pt_(0, df, mean, sd))
  }
