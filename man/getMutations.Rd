\name{getMutations}
\alias{getMutations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Mutational Data
}
\description{
Reads the mutation matrices and fasta information for each protein subunit within the quarternary structure.
}
\usage{
getMutations(mutation_files, uniprots)
}
\arguments{
  \item{mutation_files}{
A list of strings where each string is the path to a mutation matrix. A mutation matrix is a matrix of 0's (no mutation) and 1's (mutation). Each column represents a specific amino acid in the protein and each row represents an individual sample (test subject, cell line, etc). If column i in row j had a 1, that would mean that the ith amino acid for person j had a nonsynonomous mutation. As the quarternary structure can be comprised of several proteins (each with their unique uniprot id), each matrix represents the protein referenced by a specific uniprot identifier.
}
  \item{uniprots}{
A list of uniprots. The list provides the uniprot id for each of the matrices described in the \emph{mutation.data} parameter. 
}
}
\details{
The ordering in both \emph{mutation_files} and \emph{uniprots} must be the identical. For example, suppose that the quarternary structure is comprised of two proteins, A and B. If the first element in \emph{mutation_files} points to the mutation matrix for protein A, that means that the first element of \emph{uniprots}, must be a string with the uniprot id of protein A. 
}
\value{
\item{mut_tables}{A list of the mutation matrices. There should be one mutation matrix for each uniprot id in the entire assembly.}
\item{uniprots}{The uniprot ids for each of the mutation matrices. The uniprot id's are shown in the same order as the mutation matrices.}
\item{aa_counts}{The number of amino acids for each uniprot. This corresponds to the number of columns in the mutation matrix that is provided as input.}
\item{canonical_lengths}{The length of the protein as shown in the uniprot database. The uniprot ID must be available on \url{uniprot.org.}}
}

\references{
The UniProt Consortium. \emph{Activities at the Universal Protein Resource (UniProt)}. Nucleic Acids Res. 42: D191-D198 (2014).
}

\note{
To see examples of mutation matrices, please look in the /emph{extdata} folder of the package.
}

\examples{
mutation_files <- list(
					system.file("extdata","HFE_Q30201_MutationOutput.txt", package = "QuartPAC"),
          system.file("extdata","B2M_P61769_MutationOutput.txt", package = "QuartPAC")
					)
uniprots <- list("Q30201","P61769")

(mutation.data <- getMutations(mutation_files = mutation_files, uniprots = uniprots))
}
\keyword{ mutations }
