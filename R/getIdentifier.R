getIdentifier <-
function(current_identifier, ids, resSeq){
  matching_chain_rows <- which(ids$chain == current_identifier) #in case there are multiple A chains
  sequences <- lapply(1: dim(ids)[1], function(x){return (ids$seq_begin[x]: ids$seq_end[x])})
  in_list <- unlist(lapply(1:length(sequences), function(x) {ifelse(resSeq %in% sequences[[x]], x, 0)}))
  needed_row <- intersect(matching_chain_rows, in_list)
  
  if(length(needed_row) == 0)
    return(list(UNP = "NA", dbref = "NA"))
  else if (ids[needed_row,]$database == "UNP")
    return(list(UNP = ids[needed_row,]$identifier, dbref = needed_row))
  else
    return(list(UNP = "NA", dbref = "NA"))
}
