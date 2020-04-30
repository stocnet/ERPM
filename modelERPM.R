modelERGM <- function(adjacency, attributes, effects, mode, direction, startingEstimates) {
  
  
  if(mode == "cross-sectional") {

    results <- EstimationCS(adjacency, attributes, effects, direction, startingestimates, multiplicationfactor = 10, gainfactor = 0.1)
    
  } else {
    
  }
  
}