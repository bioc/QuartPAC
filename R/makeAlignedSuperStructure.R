makeAlignedSuperStructure <-
function(PDB_location, Assembly_location){
  alignment_details<- makeAlignmentMaps(PDB_location, Assembly_location) #gets a list of the alignment information, entire structure, and unprots

  #We now will make a new structure with the canonical positions appended for only the positions where canonical= pdb
  canonical_pos <- array(data = NA, dim = dim(alignment_details$structure)[1]) #set size of initial array
  
  for(i in 1:dim(alignment_details$structure)[1]){ #Iterate through array and make the code
    lookup_pos <- alignment_details$structure$absPos[i]
    alignment_row <- which(alignment_details$alignment$positional_numbering == lookup_pos)
    if(length(alignment_row) == 0) #This will occur if the no residue was aligned to the absPos. That absPo value was not able to be found in the lookup table. Thus no canonical number. 
      canonical_pos[i] <- -999
    else if(toString(alignment_details$alignment$pattern[alignment_row]) != toString(alignment_details$alignment$subject[alignment_row]))
      canonical_pos[i] <- -999 #This will be the error code -- we want this column all numeric
    else
      canonical_pos[i] <- alignment_details$alignment$canonical_numbering[alignment_row] #Match occures
  }
  
  aligned_structure <- cbind(alignment_details$structure, canonical_pos)
  setnames(aligned_structure, c(colnames(alignment_details$structure), "canonical_pos"))
  aligned_structure <- aligned_structure[which(aligned_structure$canonical_pos != -999),]  
  aligned_structure$UNP <- factor(aligned_structure$UNP)
  
  return (list(aligned_structure = aligned_structure, aa_table = alignment_details$alignment, raw_structure = alignment_details$structure, uniprots = levels(aligned_structure$UNP)))
}
