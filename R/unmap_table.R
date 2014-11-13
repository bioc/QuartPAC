unmap_table <-
function(clustering_results, positional_data){
  cluster_locations <- matrix(data = 0, nrow = dim(clustering_results)[1], ncol = 2)
  colnames(cluster_locations) <- c("start_serial", "end_serial")
  lookup_values<-clustering_results[,2:3]
  
  if(is.null(dim(lookup_values))){ #If there's only one cluster, there's no matrix so the else won't work
    cluster_locations[,1] <- positional_data$aligned_structure[lookup_values[1]]$serial
    cluster_locations[,2] <- positional_data$aligned_structure[lookup_values[2]]$serial
  }else{
    cluster_locations[,1] <- positional_data$aligned_structure[lookup_values[,1]]$serial
    cluster_locations[,2] <- positional_data$aligned_structure[lookup_values[,2]]$serial
  }
  
  clustering_results[,2:3] <- cluster_locations
  colnames(clustering_results)[1:3] <- c("AA_in_Cluster","serial_start", "serial_end")
  return(clustering_results)
}
