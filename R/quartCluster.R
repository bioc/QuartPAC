quartCluster <-
function(mutation_data, alignment, perform.ipac = "Y", perform.graphpac = "Y", perform.spacepac = "Y",
                        insertion.type = "cheapest_insertion", MultComp = "Bonferroni",
                        alpha = 0.05,
                        show.low.level.messages = "N",
                        ipac.method = "MDS", spacepac.method = "SimMax",
                        create.map = "Y", Show.Graph = "Y", Graph.Output.Path = NULL,
                        Graph.File.Name = "Map.pdf",
                        Graph.Title = "Mapping",
                        fix.start.pos = "Y",
                        numsims = 1000,
                        simMaxSpheres=3,
                        radii.vector = c(1,2,3,4,5,6,7,8,9,10),
                        OriginX = "",
                        OriginY = "",
                        OriginZ = ""){
  
  mega_mutation_matrix<- buildMegaMutationMatrix(mutation_data, alignment)
  
  #some of the legacy code GraphPAC code relies upon having this as colnames
  colnames(mega_mutation_matrix) <- unlist(lapply(1:dim(mega_mutation_matrix)[2], function(x){sprintf("V%i",x)}))
  
  custom_position_data <- cbind(
                              as.data.frame(alignment$aligned_structure$resName),
                              as.data.frame(seq(1:dim(mega_mutation_matrix)[2])),
                              as.data.frame(alignment$aligned_structure$chainID),
                              as.data.frame(alignment$aligned_structure$xCoord),
                              as.data.frame(alignment$aligned_structure$yCoord),
                              as.data.frame(alignment$aligned_structure$zCoord)
                              )
  if(OriginX == ""){
    OriginX = min(custom_position_data[, 4])
  }
  if(OriginY == ""){
    OriginY = min(custom_position_data[, 5])
  }
  if(OriginZ == ""){
    OriginZ = min(custom_position_data[, 6])
  }
  
  colnames(custom_position_data) <- c("Residue", "Can.Count","SideChain","XCoord","YCoord","ZCoord")
  
  if(perform.spacepac == "Y"){
    message("Performing SpacePAC calculations")
    if(show.low.level.messages == "Y"){
      results_spacepac <- SpaceClust(mutation.data = mega_mutation_matrix,
                                     position.matrix = custom_position_data,
                                     method = spacepac.method,
                                     numsims = numsims,
                                     simMaxSpheres = simMaxSpheres,
                                     radii.vector = radii.vector,
                                     multcomp = MultComp,
                                     alpha = alpha)
      spacepac_messages = NULL
    }else{
      spacepac_messages <- capture.output(results_spacepac <- SpaceClust(mutation.data = mega_mutation_matrix,
                                                                         position.matrix = custom_position_data,
                                                                         method = spacepac.method,
                                                                         numsims = numsims,
                                                                         simMaxSpheres = simMaxSpheres,
                                                                         radii.vector = radii.vector,
                                                                         multcomp = MultComp,
                                                                         alpha = alpha))
    }
    if(!is.null(results_spacepac$best.1.sphere)){
      results_spacepac$best.1.sphere<-unmap_sphere(1, results_spacepac$best.1.sphere, alignment)
    }
    if(!is.null(results_spacepac$best.2.sphere)){
      results_spacepac$best.2.sphere<-unmap_sphere(1:2, results_spacepac$best.2.sphere, alignment)
    }
    if(!is.null(results_spacepac$best.3.sphere)){
      results_spacepac$best.3.sphere<-unmap_sphere(1:3, results_spacepac$best.3.sphere, alignment)
    }
    
    if(results_spacepac$optimal.num.spheres == 1){
      results_spacepac$optimal.sphere <- results_spacepac$best.1.sphere
    }else if(results_spacepac$optimal.num.spheres == 2){
      results_spacepac$optimal.sphere <- results_spacepac$best.2.sphere 
    }else if(results_spacepac$optimal.num.spheres == 3){
      results_spacepac$optimal.sphere <- results_spacepac$best.3.sphere
    }else{
      results_spacepac$optimal.sphere <- NULL
    }
    
    
    
  }else{
    results_spacepac = NULL
  }
  
  
  if(perform.ipac == "Y"){
    message("Performing iPAC Calculations")
    if(show.low.level.messages == "Y"){
      results_ipac <- ClusterFind(mutation.data = mega_mutation_matrix, 
                                                                  position.data = custom_position_data,
                                                                  method = ipac.method,
                                                                  alpha = alpha,
                                                                  MultComp = MultComp,
                                                                  Include.Culled = "N",
                                                                  Include.Full = "N",
                                                                  create.map = create.map,
                                                                  Show.Graph = Show.Graph,
                                                                  Graph.Output.Path = Graph.Output.Path,
                                                                  Graph.File.Name = Graph.File.Name,
                                                                  Graph.Title = Graph.Title,
                                                                  OriginX = OriginX,
                                                                  OriginY = OriginY,
                                                                  OriginZ = OriginZ)
      ipac_messages <- NULL
    }else{
      ipac_messages <- capture.output(results_ipac <- ClusterFind(mutation.data = mega_mutation_matrix, 
                                                                  position.data = custom_position_data,
                                                                  method = ipac.method,
                                                                  alpha = alpha,
                                                                  MultComp = MultComp,
                                                                  Include.Culled = "N",
                                                                  Include.Full = "N",
                                                                  create.map = create.map,
                                                                  Show.Graph = Show.Graph,
                                                                  Graph.Output.Path = Graph.Output.Path,
                                                                  Graph.File.Name = Graph.File.Name,
                                                                  Graph.Title = Graph.Title,
                                                                  OriginX = OriginX,
                                                                  OriginY = OriginY,
                                                                  OriginZ = OriginZ))
    }
    

    if(is.null(results_ipac$Remapped)){
      results_ipac = NULL
    }else{
      results_ipac <- unmap_table(results_ipac$Remapped, alignment)
    }
  }else{
    results_ipac = NULL
  }
  if(perform.graphpac == "Y"){
    message("Performing GraphPAC Calculations")
    if(show.low.level.messages == "Y"){
      results_graphpac<- GraphClust(mutation.data = mega_mutation_matrix,
                                    position.data = custom_position_data,
                                    insertion.type = insertion.type, 
                                    alpha = alpha, 
                                    MultComp = MultComp,
                                    fix.start.pos = fix.start.pos,
                                    Include.Culled = "N",
                                    Include.Full = "N")
      graphpac_messages <- NULL
    }else{
      graphpac_messages <- capture.output(results_graphpac<- GraphClust(mutation.data = mega_mutation_matrix,
                                                                        position.data = custom_position_data,
                                                                        insertion.type = insertion.type, 
                                                                        alpha = alpha, 
                                                                        MultComp = MultComp,
                                                                        fix.start.pos = fix.start.pos,
                                                                        Include.Culled = "N",
                                                                        Include.Full = "N"))
    }

    if(is.null(results_graphpac$Remapped)){
      results_graphpac = NULL
    }else{
      results_graphpac$Remapped <- unmap_table(results_graphpac$Remapped, alignment)
      results_graphpac$candidate.path <- unmap_sequence(results_graphpac$candidate.path, alignment)
      results_graphpac <- list(clusters = results_graphpac$Remapped, candidate.path = results_graphpac$candidate.path, path.distance = results_graphpac$path.distance, linear.path.distance = results_graphpac$linear.path.distance)
    }
    
  }else{
    results_graphpac = NULL
  }

  return(list(ipac = results_ipac, graphpac = results_graphpac, spacepac = results_spacepac, ipac_messages = ipac_messages,
              graphpac_messages = graphpac_messages, spacepac_messages=spacepac_messages))
}
