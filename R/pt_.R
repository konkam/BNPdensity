pt_ <-
  function(x, df, mean, sd) {
    pt((x - mean) / sd, df, ncp = 0)
  }
