alignPositions <-
function(canonical.sequence, position.df, atom = "CA",  
                           patternQuality = PhredQuality(22L), subjectQuality = PhredQuality(22L),
                           type = "global-local", substitutionMatrix =NULL, fuzzyMatrix = NULL,
                           gapOpening = -10, gapExtension = -4, scoreOnly = FALSE){
  
  canonical.sequence <- paste(canonical.sequence, collapse="")
  position.sequence <- paste(unlist(lapply(position.df$resName, get.SingleLetterCode)), collapse="")
  
  alignment <- pairwiseAlignment(pattern = canonical.sequence, subject = position.sequence, patternQuality=patternQuality,
                                 subjectQuality=subjectQuality,type = type, substitutionMatrix= substitutionMatrix,
                                 fuzzyMatrix=fuzzyMatrix,gapOpening=gapOpening,gapExtension=gapExtension,
                                 scoreOnly=scoreOnly)
  return(alignment)
}
