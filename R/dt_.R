dt_ <-
  function(x, df, mean, sd) {
    dt((x - mean) / sd, df, ncp = 0) / sd
  }
