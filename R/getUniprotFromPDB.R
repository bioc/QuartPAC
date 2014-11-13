getUniprotFromPDB <-
function(data, defs){
    dbref_rows <- which(unlist(lapply(data, function(x){substring(x,defs$recordName$start,defs$recordName$end)})) == "DBREF ") #Note the extra space
    identifiers<- lapply(data[dbref_rows], function(x){
    chain <- trim.both(substring(x, defs$chainID$start, defs$chainID$end))
    database <- trim.both(substring(x, defs$database$start, defs$database$end))
    identifier <- trim.both(substring(x, defs$dbAccession$start, defs$dbAccession$end))
    seq_begin <- as.numeric(trim.both(substring(x, defs$seqBegin$start, defs$seqBegin$end)))
    seq_end <- as.numeric(trim.both(substring(x, defs$seqEnd$start, defs$seqEnd$end)))
    return( cbind(data.frame(database), data.frame(chain), data.frame(identifier), data.frame(seq_begin), data.frame(seq_end)))
    
    # data.frame(cbind (database, chain, identifier, seq_begin, seq_end)))
    
  })
  identifiers <- rbindlist(identifiers)
  return(identifiers)
}
