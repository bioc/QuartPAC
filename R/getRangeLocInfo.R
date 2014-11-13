getRangeLocInfo <-
function(serial_start, serial_end, positional_data){
  uniprot_list <- lapply(positional_data$uniprots, function(x){list(serial_start = NA, serial_end=NA)})
  names(uniprot_list) <- positional_data$uniprots
  rows_in_range <- which(positional_data$aligned_structure$serial == serial_start):which(positional_data$aligned_structure$serial == serial_end) 
  serials_in_range <- positional_data$aligned_structure$serial[rows_in_range]
  uniprots_in_range <- positional_data$aligned_structure$UNP[rows_in_range]
  chainIDs_in_range <- positional_data$aligned_structure$chainID[rows_in_range]
  resSeqs_in_range <- positional_data$aligned_structure$resSeq[rows_in_range]
  long_mapped_data <- cbind(as.data.frame(uniprots_in_range), as.data.frame(chainIDs_in_range), as.data.frame(resSeqs_in_range), as.data.frame(serials_in_range))
  colnames(long_mapped_data) <- c("Uniprot", "Chain", "resSeq", "Serial")
  return(mapped_range = do.call(data.frame, aggregate(cbind(resSeq, Serial) ~ Uniprots + Chain, data = long_mapped_data, function(x) c(start = min(x), end = max(x)))))
}
