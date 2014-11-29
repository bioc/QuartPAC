\name{makeAlignedSuperStructure}
\alias{makeAlignedSuperStructure}
\title{
Create Protein Assembly Information
}
\description{
Reads the information in the PDB files to build the quaternary structure.
}
\usage{
makeAlignedSuperStructure(PDB_location, Assembly_location)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PDB_location}{
A string pointing to the location of the *.pdb file. The pdb file provides uniprot and structural information for each residue.
}
  \item{Assembly_location}{
A string pointing to the lcoation of the *.pdb1 file. The pdb1 file provides the positional information of all the residues in the structure when they are in the asembly.
}
}
\value{

\item{aligned_structure}{The aligned structure. For each residue, you get the XYZ coordinate as well as other structural information. The absPos column in the \emph{aligned_structure} variable matches the absPos column in the \emph{raw_structure} and is used for debugging purposes. The \emph{canonical_pos} column represents the residue number in the fasta sequence for the subunit with the specified uniprot id. The \emph{canonical_pos} column can repeat. For instance, if there are two proteins, A and B, in the assembly, both proteins may have residue #5 in them. The \emph{absPos} column will not repeat.}
\item{aa_table}{A table showing the pairwise alignment between the residue sequence as specified by the pdb file and the residue sequence specified by the fasta sequence in the uniprot database. This table is mainly used for debugging.}
\item{raw_structure}{The raw structure as read in from the *.pdb1 file after some cleaning. For instance, row's with NA for uniprot values are dropped. Further, only the carbon-alpha rows are kept.}
}
\note{
This method reads fasta information from \url{www.uniprot.org}. Thus an internet connection is required in order to execute properly.

This functionality is still in beta. The user is encouraged to check the alignment created by the method. 
}

\examples{
#read the pdb file
pdb.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb"
assembly.location <- "http://www.rcsb.org/pdb/files/1A6Z.pdb1"


(alignment.results <- makeAlignedSuperStructure(pdb.location, assembly.location))
}
\keyword{ structure }
