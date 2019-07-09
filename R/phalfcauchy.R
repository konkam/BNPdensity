phalfcauchy <-
  function(q, location = 0, scale = 1) {
    ifelse(q < 0, 0, 1) * (pcauchy(q, location, scale) - pcauchy(
      0,
      location, scale
    )) / (1 - pcauchy(0, location, scale))
  }
