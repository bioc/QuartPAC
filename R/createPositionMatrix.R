createPositionMatrix <-
function(PDB_location, Assembly_location, requiredAtom = "CA"){
  pdb_file <- readLines(PDB_location)
  assembly_file <- readLines(Assembly_location)
  
  file.data <- loadDefs()
  PDB_defs = file.data$PDB_defs
  DBREF_defs = file.data$DBREF_defs
  
  #find cut positions in structure
  model_cuts<- findModelCuts(assembly_file)
  ##get chain identifiers
  ids <- getUniprotFromPDB(pdb_file, DBREF_defs)
  
  structure <- list()
  ##Create one structure per model
  structure <- lapply(1:model_cuts$num_protomers, function(x){extractAtoms(assembly_file[model_cuts$model_starts[x]: model_cuts$model_ends[x]], PDB_defs, ids)} )
  structure <- lapply(1:model_cuts$num_protomers, function(protomer){cbind(structure[[protomer]], protomer)}) #attach protomer numbers location
  
  structure<- rbindlist(structure) #make superstructure
  structure$dbref <- as.factor(structure$dbref) #convert this to a factor - this fixes a bug when no "NA's" in dbref's (perfect lookup to dbref table). Later on, makeAlignmentMaps.R expects a factor.
  
  structure <- structure[which(structure$atom == requiredAtom),] #select only the desired atoms
  
  #cut out rows with bad uniprot
  bad_rows <- which(structure$UNP == "NA")
  if(length(bad_rows) != 0)
    structure <- structure[-bad_rows, ]
  
  #attach raw number for positional representation
  structure <- cbind(structure, absPos = 1:dim(structure)[1]) #Attach RawNumber
  
  return (structure)
  
}
