test_quartCluster <- function() {
  
  ##############
  #check the mutational data is outputted as expected for typical example
  #############
  mutation_files <- list(
    system.file("extdata","HFE_Q30201_MutationOutput.txt", package = "QuartPAC"),
    system.file("extdata","B2M_P61769_MutationOutput.txt", package = "QuartPAC")
  )
  uniprots <- list("Q30201","P61769")
  mutation.data <- getMutations(mutation_files = mutation_files, uniprots = uniprots)
 
  #Begin checks
  checkEquals(length(mutation.data), 6)
  checkEquals(dim(mutation.data$mut_tables[[1]])[1], 4)
  checkEquals(dim(mutation.data$mut_tables[[1]])[2], 348)
  checkEquals(dim(mutation.data$mut_tables[[2]])[1], 7)
  checkEquals(dim(mutation.data$mut_tables[[2]])[2], 119)
  checkEquals(sum(mutation.data$mut_tables[[1]]), 4)
  checkEquals(sum(mutation.data$mut_tables[[2]]), 7)
  checkEquals(mutation.data$canonical_lengths[1], 348)
  checkEquals(mutation.data$canonical_lengths[2], 119)
  checkEquals(mutation.data$aa_counts[1], 348)
  checkEquals(mutation.data$aa_counts[2], 119)
  checkEquals(mutation.data$uniprots[1], "Q30201")
  checkEquals(mutation.data$uniprots[2], "P61769")
  
  
  #########
  #check the structural data data is outputted as expected for typical example
  #########
  pdb.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb"
  assembly.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb1"
  structural.data <- makeAlignedSuperStructure(pdb.location, assembly.location)
  
  #Begin Checks
  checkEquals(length(structural.data), 4)
  checkEquals(dim(structural.data$aligned_structure)[1], 371)
  checkEquals(dim(structural.data$aligned_structure)[2], 20)
  checkEquals(dim(structural.data$aligned_structure)[1], 371)
  checkEquals(dim(structural.data$aligned_structure)[2], 20)
  checkEquals(structural.data$uniprots[1], "Q30201")
  checkEquals(structural.data$uniprots[2], "P61769")
  #check uniprots match in mutations and structural data
  checkEquals(mutation.data$uniprots[1], structural.data$uniprots[1])
  checkEquals(mutation.data$uniprots[2], structural.data$uniprots[2])
  
  
  ####
  #Run the algorithm and check the results
  ####
  quart_results <- quartCluster(mutation.data, structural.data, perform.ipac = "Y", perform.graphpac = "Y",
                                perform.spacepac = "Y", create.map = "N", MultComp = "None",
                                alpha = 1, radii.vector = c(1:2), show.low.level.messages = "Y")
  
  #BeginChecks
  checkEquals(length(quart_results), 6)
  checkEquals(dim(quart_results$ipac)[1], 35)
  checkEquals(dim(quart_results$ipac)[2], 5)
  checkEquals(dim(quart_results$graphpac$clusters)[1], 35)
  checkEquals(dim(quart_results$graphpac$clusters)[2], 5)
  checkEquals(length(quart_results$graphpac$candidate.path), 371)
  checkTrue(quart_results$spacepac$optimal.num.spheres>0)
  checkTrue(quart_results$spacepac$optimal.radius > 0)
  checkTrue(quart_results$spacepac$p.value <= 1)
  checkTrue(quart_results$graphpac$path.distance > 0)
  checkTrue(is.numeric(quart_results$graphpac$path.distance))
  checkTrue(quart_results$graphpac$linear.path.distance > 0)
  checkTrue(is.numeric(quart_results$graphpac$linear.path.distance))
  checkTrue(is.null(quart_results$ipac_messages))
  checkTrue(is.null(quart_results$graphpac_messages))
  checkTrue(is.null(quart_results$spacepac_messages))
}