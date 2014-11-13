extractAtoms <-
function (data, cuts, ids){
  atom_rows <- which(unlist(lapply(data, function(x){substring(x,cuts$recordName$start,cuts$recordName$end)})) == "ATOM  ")
  positions <- rbindlist(lapply(data[atom_rows], extractAtomRow, cuts, ids))
  return(positions)
}
