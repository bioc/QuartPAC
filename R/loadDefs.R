loadDefs <-
function(){
  PDB_defs = list(
    recordName = list(start =1, end = 6), 
    serial = list(start = 7, end = 11),
    atom = list(start = 13, end =16),
    altLoc = list(start = 17, end = 17),
    resName = list(start = 18, end= 20),
    chainID = list(start = 22, end = 22),
    resSeq = list(start = 23, end = 26),
    iCode = list(start = 27, end = 27),
    x = list(start= 31, end = 38),
    y = list(start = 39, end = 46),
    z = list(start = 47, end = 54),
    occupancy = list(start = 55, end= 60), 
    tempFactor = list(start = 61, end = 66),
    element = list(start = 77, end = 78),
    charge = list(start = 79, end = 80))
  
  DBREF_defs = list(
    recordName = list(start =1, end =6),
    idCode = list(start = 8, end = 11),
    chainID = list(start = 13, end = 13),
    seqBegin = list(start = 15, end = 18),
    insertBegin = list(start = 19, end = 19),
    seqEnd = list(start = 21, end = 24),
    insertEnd = list(start = 25, end = 25),
    database = list(start = 27, end = 32),
    dbAccession = list(start= 34, end = 41),
    dbIdCode = list(start = 43, end = 54),
    dbseqBegin = list(start = 56, end = 60),
    idbnsBegin = list(start = 61, end =61),
    dbseqEnd = list(start = 63, end = 67),
    dbinsEnd = list(start = 68, end = 68))
  
  return (list(PDB_defs = PDB_defs, DBREF_defs=DBREF_defs))
  
}
