\name{ParallelStructure-package}
\alias{ParallelStructure-package}
\alias{ParallelStructure}
\docType{package}
\title{
R FRAMEWORK TO RUN PARALLEL JOBS IN  STRUCTURE
}
\description{
Contains 

 - \link{MPI_structure} and   \link{parallel_structure}: two R functions that calls STRUCTURE automatically for a set of predefined jobs. Jobs are distributed among all available computing units (cores or cpu) in order to make the best use of multi-core computer while running analysis of large data sets in STRUCTURE. For the user perspective  \link{MPI_structure} and   \link{parallel_structure} are equivalent and work with the same set of input. Performance might vary slightly as they rely on two different method for distributing the jobs among multiple processors. 
 
  \link{MPI_structure} relies on "Rmpi" package which is available from the CRAN package repository for MacOS but requires manual installation for Windows version of R.
  
   \link{parallel_structure}: Relies on R  \link{parallel-package} which is available for both Mac and Windows version of R since version 2.14. However the package is under development and have some limitation (should not be used in GUI  or embedded environments).
 
}
\details{
\tabular{ll}{
Package: \tab ParallelStructure\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-07\cr
License: \tab free\cr
LazyLoad: \tab yes\cr
}

  \link{MPI_structure}  uses Rmpi package to pilot serial STRUCTURE jobs in parallel  
MPI_structure requires:

1- the R package Rmpi : http://www.stats.uwo.ca/faculty/yu/Rmpi/

	-MacOS: A functional version of Rmpi is available from the package repository CRAN
		(otherwise see: http://www.stats.uwo.ca/faculty/yu/Rmpi/)
		
	- Windows and LINUX manual installation is required for RMPI, see specific instruction on http://www.stats.uwo.ca/faculty/yu/Rmpi/

2- the software STRUCTURE installed: http://pritch.bsd.uchicago.edu/structure.html  
}
\author{
Francois Besnier, Kevin Glover

francois.besnier@imr.no

}
\references{

Besnier F, Glover KA ParallelStructure: a R package to distribute parallel runs of the population genetics program STRUCTURE on multi-core computers. To be submitted

RMPI website: http://www.stats.uwo.ca/faculty/yu/Rmpi/

\link{parallel-package} 
 
 STRUCTURE website:  http://pritch.bsd.uchicago.edu/software.html
}
