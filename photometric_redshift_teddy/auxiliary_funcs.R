normalize_density <- function(bin_size, estimates,
                              tolerance = 1e-3, max_iters = 500) {
  area <- bin_size * sum(pmax(estimates, 0.0))
  if (area == 0) {
    estimates[] <- 1 / (length(estimates) * bin_size)
    return(estimates)
  }
  else if (area < 1) {
    estimates[estimates > 0] <- estimates[estimates > 0] / area
    return(estimates)
  }
  upper <- max(estimates)
  lower <- 0.0
  middle <- (upper + lower) / 2
  
  iter <- 1
  while (iter <= max_iters) {
    iter <- iter + 1
    
    density <- pmax(estimates - middle, 0.0)
    area <- bin_size * sum(density)
    
    if (abs(area - 1) < tolerance) {
      break
    }
    
    if (area > 1) {
      lower <- middle
    } else {
      upper <- middle
    }
    middle <- (upper + lower) / 2
  }
  
  return(pmax(estimates - middle, 0.0))
}
