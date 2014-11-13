getSequenceLocInfo <-
function(serial_vector, positional_data){
  uniprot_list <- lapply(positional_data$uniprots, function(x){list(serial_start = NA, serial_end=NA)})
  names(uniprot_list) <- positional_data$uniprots
  rows_in_range <- unlist(lapply(serial_vector, function(x){which(positional_data$aligned_structure$serial ==x)}))
  serials_in_range <- positional_data$aligned_structure$serial[rows_in_range]
  uniprots_in_range <- positional_data$aligned_structure$UNP[rows_in_range]
  chainIDs_in_range <- positional_data$aligned_structure$chainID[rows_in_range]
  resSeqs_in_range <- positional_data$aligned_structure$resSeq[rows_in_range]
  long_mapped_data <- cbind(as.data.frame(uniprots_in_range), as.data.frame(chainIDs_in_range), as.data.frame(resSeqs_in_range), as.data.frame(serials_in_range))
  colnames(long_mapped_data) <- c("Uniprot", "Chain", "resSeq", "Serial")
  return(long_mapped_data)
  
}
