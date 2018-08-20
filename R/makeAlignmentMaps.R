makeAlignmentMaps <-
function(PDB_location, Assembly_location){
  protein_structure <- createPositionMatrix(PDB_location, Assembly_location)
  #protein_structure$UNP <- factor(protein_structure$UNP) # in case of misread.
  dbrefs <- levels(protein_structure$dbref)
  dbrefs <- as.numeric(dbrefs[which(dbrefs != "NA")])
  uniprots <- unlist(lapply(dbrefs, function(x){protein_structure[which(protein_structure$dbref ==x),]$UNP[1]}))
  alignment_maps <- lapply(dbrefs, function(x){ 
                        current_fasta <- sprintf("https://www.uniprot.org/uniprot/%s.fasta",toString(uniprots[x]) )
                        current_meta_structure <- protein_structure[which(protein_structure$dbref == dbrefs[x]),]
    
                        #Now we break up the current_structure by protomer -- in case one dbref has multiple protomers
                        structures <- lapply(unique(current_meta_structure$protomer), function(y){
                                              return (current_meta_structure[which(current_meta_structure$protomer == y),])
                                      })
                        protomer_alignment_maps<- rbindlist(lapply(structures, constructCanonicalNumbering, current_fasta))
                        
    return (protomer_alignment_maps)
  })
  return (list(alignment = rbindlist(alignment_maps), structure = protein_structure, uniprots = uniprots))
}
