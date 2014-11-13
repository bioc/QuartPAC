constructCanonicalNumbering <-
function(structure, fasta){
  current_seq <- get.AASeq(fasta)
  alignment <- alignPositions(current_seq, structure)
  
  alignment_matrix <- cbind(
    as.data.frame(unlist(strsplit(toString(pattern(alignment)),""))),
    as.data.frame(unlist(strsplit(toString(subject(alignment)),""))))
  colnames(alignment_matrix) <- c("pattern", "subject")
  
  canonical_numbering <- positional_numbering <- array(data = NA, length(dim(alignment_matrix)[1]))
  
  start_canonical <- start(pattern(alignment))
  start_positional <- start(subject(alignment))
  current_count_canonical <- 0
  current_count_positional <- structure$absPos[1]-1
  
  
  for(i in 1: dim(alignment_matrix)[1]){
    if(alignment_matrix[i,1]!="-"){
      canonical_numbering[i] <- start_canonical+current_count_canonical
      current_count_canonical<- current_count_canonical+1
    }
    if(alignment_matrix[i,2]!="-"){
      positional_numbering[i] <- start_positional+current_count_positional
      current_count_positional<- current_count_positional+1
    }
  }
  alignment_matrix <- cbind(alignment_matrix, canonical_numbering, positional_numbering)
  return (alignment_matrix)
}
