unmap_row <-
function(row, sphere_numbers, positional_data){
  
  for(sphere_number in sphere_numbers){
  
    if(length(row$Center1) == 0){ ##must be in best.1.sphere
      sphere_number = "" #best.1.sphere just has Center/Start/End -- thus paste should add nothing.
    }
    
    if(sphere_number ==1){
      if(length(unlist(row$Intersection)) !=0){
        row$Intersection  <- "Error! Should be 0."
      }
    }
  
    center_serials <- unmap_sequence(unlist(row[[paste("Center",sphere_number,sep="")]]), positional_data)
    row[[paste("Center",sphere_number,sep="")]] <- if(length(center_serials)==0){NULL}else{list(center_serials)}
  
    start_serials <-unmap_sequence(unlist(row[[paste("Start",sphere_number,sep="")]]), positional_data)
    row[[paste("Start",sphere_number,sep="")]] <- if(length(start_serials)==0){NULL}else{list(start_serials)}
  
    end_serials <- unmap_sequence(unlist(row[[paste("End",sphere_number,sep="")]]), positional_data)
    row[[paste("End",sphere_number,sep="")]] <- if(length(end_serials) ==0){NULL}else{list(end_serials)}
  
    position_serials <- unmap_sequence(unlist(row[[paste("Positions",sphere_number,sep="")]]), positional_data)
    row[[paste("Positions",sphere_number,sep="")]] <- if(length(position_serials) ==0){NULL}else{list(position_serials)}
  
  
    within.range_serials <- unmap_sequence(unlist(row[[paste("Within.Range",sphere_number,sep="")]]), positional_data)
    row[[paste("Within.Range",sphere_number,sep="")]] <- if(length(within.range_serials)==0){NULL}else{list(within.range_serials)}
  }
  
  return(row)
}
