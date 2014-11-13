buildMegaMutationMatrix <-
function(mutation_data, alignment){
  #check lengths
  if(sum(mutation_data$aa_counts != mutation_data$canonical_lengths)!=0){
    stop("Error: Your mutation matrices do not have the expeted number of columns. For each Uniprot, the number of columns of the matrix must match the length of the protein.")
  }
  
  #check number of mutations
  if(sum(mutation_data$mutation_counts) < 2){
    stop("Error: You have only 1 mutation for this structure! No clustering possible.")
  }
  
  #check that all the uniprots are there.
  if(length(setdiff(mutation_data$uniprots, alignment$uniprots)) != 0){
    missing_in_mutations<- toString(setdiff(alignment$uniprots, mutation_data$uniprots))
    missing_in_alignment <- toString(setdiff(mutation_data$uniprots, alignment$uniprots))
    error_message <- sprintf("Error: Not all uniprots are in both the mutational and positional data. \n 
                            Uniprots in alignment matrices but not in mutation matrices: %s \n
                            Uniprots in mutation matrices but not in alignment matrices: %s", missing_in_mutations, missing_in_alignment)
    stop(error_message)
  }
  
  #Construct sequential matrix that will later be rearranged
  big_matrix <- Reduce(matrix_combine, mutation_data$mut_tables)
  
  
  #Now we construct the correct ordered matrix
  big_muts <- matrix(data= NA, nrow = dim(big_matrix)[1], ncol = dim(alignment$aligned_structure)[1])
  
  
  #construct column headers
  for(i in 1:length(mutation_data$aa_counts)){
    if(i == 1){
      column_headers <- list( list(start = 1, end = mutation_data$aa_counts[i]))
    }else{
      column_headers[[length(column_headers)+1]] <- list( start =column_headers[[i-1]]$end+1, end =column_headers[[i-1]]$end+1+mutation_data$aa_counts[i]-1)
    }
  }
  
  
  for(i in 1:length(column_headers)){
    current_uniprot<- mutation_data$uniprots[i]
    colnames(big_matrix)[column_headers[[i]]$start:column_headers[[i]]$end] <- current_uniprot
  }
  
  #rearrange_matrix
  for(i in 1:(dim(big_muts)[2])){
    current_unp <- alignment$aligned_structure$UNP[i]
    required_position <- alignment$aligned_structure$canonical_pos[i]
    big_muts[,i]<-big_matrix[,which(colnames(big_matrix)==current_unp)[required_position]]
  }
  
  return(big_muts)
}
