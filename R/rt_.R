rt_ <-
  function(n, df, mean, sd) {
    mean + sd * rt(n, df, ncp = 0)
  }
