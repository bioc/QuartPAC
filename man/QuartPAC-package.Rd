\name{QuartPAC-package}
\alias{QuartPAC-package}
\alias{QuartPAC}
\docType{package}
\title{
Identifying mutational clusters while incorporating protein quaternary structure.
}
\description{
\pkg{QuartPAC} is a companion package to \pkg{iPAC}, \pkg{GraphPAC} and \pkg{iPAC}. It allows one to use the methodologies proposed by each of those packages to be applied to the protein quaternary structure.
}
\details{
QuartPAC is designed to identify mutational clustering in 3D space when looking at the entire quaternary protein structure. It does this by applying the algorithms proposed in \pkg{iPAC}, \pkg{GraphPAC} and \pkg{SpacePAC} over the entire assembly after correctly preprocessing the mutational and positional data. QuartPAC typically follows a three step process. Step 1 involves reading the mutational data - see \code{getMutations} for more information. Step 2 involves reading in the structural information to create the protein assembly and is explained in \code{makeAlignedSuperStructure}. Finally, Step 3 performs the statistical analysis and reports the significant p-values -- see \code{QuartCluster} for more information.

The clustering results give the \emph{serial number} values from the *.pdb1 file. 
}
\author{
Gregory Ryslik, Yuwei Cheng, Hongyu Zhao.
Maintainer: Gregory Ryslik <gregory.ryslik@yale.edu>
}

\references{
Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein Amino acid Clustering. R package version 1.8.0. http://www.bioconductor.org/.

Gregory Ryslik and Hongyu Zhao (2013). GraphPAC: Identification of Mutational Clusters in Proteins via a Graph Theoretical Approach. R Package version 1.6.0 http://www.bioconductor.org/.

Gregory Ryslik and Hongyu Zhao (2013). SpacePAC: Identification of Mutational Clusters in 3D Protein Space via Simulation. R package version 1.2.0. http://www.bioconductor.org/.

The UniProt Consortium. Activities at the Universal Protein Resource (UniProt). Nucleic Acids Res. 42: D191-D198 (2014).
}
\keyword{protein}
\keyword{quaternary}
\keyword{clustering}
\examples{

#read the mutational data
mutation_files <- list(
        system.file("extdata","HFE_Q30201_MutationOutput.txt", package = "QuartPAC"),
        system.file("extdata","B2M_P61769_MutationOutput.txt", package = "QuartPAC")
					)
uniprots <- list("Q30201","P61769")
mutation.data <- getMutations(mutation_files = mutation_files, uniprots = uniprots)

#read the pdb file
pdb.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb"
assembly.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb1"
structural.data <- makeAlignedSuperStructure(pdb.location, assembly.location)

#Perform Analysis
#We use a very high alpha level here with no multiple comparison adjustment
#to make sure that each method provides shows a result.
#Lower alpha cut offs are typically used.
(quart_results <- quartCluster(mutation.data, structural.data, perform.ipac = "Y", perform.graphpac = "Y",
                              perform.spacepac = "Y", create.map = "N", MultComp = "None",
                              alpha = .3, radii.vector = c(1:3), show.low.level.messages = "Y"))
}
