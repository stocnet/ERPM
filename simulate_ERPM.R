simulate_ERPM <- function(nodes, 
                          objects, 
                          effects, 
                          minestimates, 
                          maxestimates,
                          burnin = 30, 
                          thining = 10,
                          num.steps = 100,
                          neighborhood = 1) {
  
  num.nodes <- nrow(nodes)
  first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
  first.partition <- order_groupids(first.partition)
  
  # if only one chain to draw
  if(length(maxestimates) == 0) {
    
    chain <- draw_Metropolis(startingestimates, first.partition, nodes, effects, objects, burnin, thining, num.steps, "normalized", neighborhood)
    draws <- chain$draws
    
    df <- data.frame(x,y1,y2)
    
    ggplot(df, aes(x)) +                    # basic graphical object
      geom_point(aes(y=y1), colour="red") +  # first layer
      geom_point(aes(y=y2), colour="green")  # second layer
    
    
  # else plot the   
  } else {
    
  }
 
  
  return
  
}