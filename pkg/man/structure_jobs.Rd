\name{structure_jobs}
\alias{structure_jobs}
\docType{data}
\title{
  JOB LIST FOR STRUCTURE RUN
}
\description{
Each line corresponds to an individual job. Wile STRUCTURE input format requires a different dataset for each set of population, \code{ParallelStructure} offers the possibility for the user to work from a large input file containing all the populations one might need to analyze in STRUCTURE. Then, a set of "jobs" are defined in which all or only a subset of population can be included. 
This avoid making a different input file for each population subset.

For each job, the user defines the set of populations to be included, STRUCTURE parameter K, burin and number of iteration. If all population in the data must be analyzed pairwise (all versus all), the list of populations for the given job ca be replaced by "pairwaise.matrix" (see job T11 in example joblist)

}
\usage{data(structure_jobs)}
\format{
 The format is:
 One line per job with following columns:
 1- job id
 
 2-list of population to be included (coma separated). If the list of population is replaced by "pairwise.matrix", an independent STRUCTURE run will be performed for each pair of population in the dataset
 
 3- K: number of cluster in the data
 
 4- burn-in value
 
 5- number of iterations
 
}
\examples{
data(structure_jobs)
}
\keyword{datasets}
