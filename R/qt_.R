qt_ <-
  function(p, df, mean, sd) {
    sd * qt(p, df, ncp = 0) + mean
  }
