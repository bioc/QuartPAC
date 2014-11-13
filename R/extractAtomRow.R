extractAtomRow <-
function(row, cuts, ids){
  recordName <- trim.both(substring(row, cuts$recordName$start, cuts$recordName$end))
  serial <- as.integer(trim.both(substring(row, cuts$serial$start, cuts$serial$end)))
  atom <- trim.both(substring(row, cuts$atom$start, cuts$atom$end))
  altLoc <- trim.both(substring(row, cuts$altLoc$start, cuts$altLoc$end))
  resName <- trim.both(substring(row, cuts$resName$start, cuts$resName$end))
  chainID <- trim.both(substring(row, cuts$chainID$start, cuts$chainID$end))
  resSeq <- as.integer(trim.both(substring(row, cuts$resSeq$start, cuts$resSeq$end)))
  iCode <- trim.both(substring(row, cuts$iCode$start, cuts$iCode$end))
  xCoord <- as.numeric(trim.both(substring(row, cuts$x$start, cuts$x$end)))
  yCoord <- as.numeric(trim.both(substring(row, cuts$y$start, cuts$y$end)))
  zCoord <- as.numeric(trim.both(substring(row, cuts$z$start, cuts$z$end)))
  occupancy <- as.numeric(trim.both(substring(row, cuts$occupancy$start, cuts$occupancy$end)))
  tempFactor <- as.numeric(trim.both(substring(row, cuts$tempFactor$start, cuts$tempFactor$end)))
  element <- trim.both(substring(row, cuts$element$start, cuts$element$end))
  charge <- trim.both(substring(row, cuts$charge$start, cuts$charge$end))
  
  identifier <- getIdentifier(chainID, ids, resSeq)
  UNP <- identifier$UNP
  dbref <- identifier$dbref
  row_result <- data.frame(recordName, serial, atom, altLoc, resName, chainID, resSeq, iCode, xCoord, yCoord, zCoord, occupancy, tempFactor, element, charge, UNP, dbref)
  
  return(row_result)
}
