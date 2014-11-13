findModelCuts <-
function(data){
  model_starts <- which(unlist(lapply(data, function(x){substring(x, 1, 5) })) == "MODEL")
  model_ends <- which(unlist(lapply(data, function(x){substring(x,1,6)})) == "ENDMDL")
  return(list(model_starts = model_starts+1, model_ends = model_ends-1, num_protomers = length(model_starts)))  #+1/-1 so that we get the actual model positions
}
