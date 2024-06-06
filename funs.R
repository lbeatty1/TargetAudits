grid_emis_objective = function(alphas){
  
  #merge in probs (alphas) with emissions
  prob_by_bin[,prob:=alphas]
  
  grid_probs = emis_grid[prob_by_bin, on=.(count_bin, emis_bin)]
  
  well_prob = emis[grid_probs, on=.(ID)]
  
  #return expected discovered emissions
  well_prob[,detected:=ER*prob]

  detected_emis = sum(well_prob$detected)
  print(paste("detected emissions:", detected_emis))
  return(-detected_emis) #optim minimizes
}



grid_emis_constraint = function(alphas){
  
  #merge in probs (alphas) with emissions
  prob_by_bin[,prob:=alphas]
  
  grid_probs = emis_grid[prob_by_bin, on=.(count_bin, emis_bin)]
  
  well_prob = emis[grid_probs, on=.(ID)]
  ninspections= sum(well_prob$prob)
  print(paste('number of inspections:', ninspections))
  return(ninspections-budget)
  
}
