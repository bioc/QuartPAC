getMutations <-
function(mutation_files, uniprots){
  mutation_data <- lapply(mutation_files, function(x){read.table(x, header= FALSE, sep = " ")})
  mutation_data <- lapply(mutation_data, function(x){as.matrix(x)}) #convert into matrices
  fastas <- lapply(uniprots, function(x){ get.AASeq(sprintf("http://www.uniprot.org/uniprot/%s.fasta",x))})
  aa_counts <- unlist(lapply(mutation_data, function(x){dim(x)[2]}))
  mutation_counts <- unlist(lapply(mutation_data, function(x){sum(x)}))
  canonical_lengths <- unlist(lapply(fastas, function(x){length(x)}))
  result <- list(mut_tables = mutation_data, fastas = fastas, uniprots = unlist(uniprots), 
                 aa_counts = aa_counts, mutation_counts = mutation_counts, canonical_lengths =canonical_lengths )
  return(result)
}
