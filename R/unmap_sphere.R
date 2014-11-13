unmap_sphere <-
function(sphere_numbers, sphere_data, positional_data){
  result <- lapply(1:dim(sphere_data)[1], function(x){unmap_row(sphere_data[x,],sphere_numbers,positional_data)})
  result <- do.call("rbind", result)
  return(result)
}
