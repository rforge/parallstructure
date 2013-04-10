\name{MPI_structure}
\alias{MPI_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FUNCTION TO RUN PARALLEL JOBS IN  STRUCTURE
}
\description{
 \code{MPI_structure} is a R function that calls STRUCTURE automatically for a set of predefined jobs. Jobs are distributed among all available computing units (cores or cpu) in order to make the best use of multi-core computer while running analysis of large data sets in STRUCTURE. Distribution of jobs to multiple cpu relies on the \bold{Rmpi} package.

TWO inputs are give to the MPI_structure function:

1-The function argument in R

2-The job specifications written in a "joblist" file, see example joblist \link{structure_jobs}.


}
\usage{
MPI_structure(joblist = NULL, n_cpu = NULL, structure_path = Mac_path, infile = NULL, outpath = NULL,

 numinds = NULL, numloci = NULL, plot_output = 1, label = 1, popdata = 1, popflag = 0, locdata = 0, 
 
 phenotypes = 0, markernames = 0, mapdist = 0, onerowperind = 0, phaseinfo = 0, recessivealleles = 0, 
 
 phased = 0, extracol = 0, missing = -9, ploidy = 2, noadmix = 0, linkage = 0, usepopinfo = 0, locprior = 0,
 
  inferalpha = 1, alpha = 1, popalphas = 0, unifprioralpha = 1, alphamax = 10, alphapropsd = 0.025, freqscorr = 1, 
  
  onefst = 0, fpriormean = 0.01, fpriorsd = 0.05, inferlambda = 0, lambda = 1, computeprob = 1, 
  
  pfromflagonly = 0, ancestdist = 0, startatpopinfo = 0, metrofreq = 10, updatefreq = 1, printqhat = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{joblist}{
  Name of the file where list of jobs is stored, see example data: joblist1
}
  \item{n_cpu}{
   number of cpu cores to be used
}
  \item{structure_path}{
    Location of the executable command line STRUCTURE program.  "/Applications/Structure.app/Contents/Resources/Java/bin/" for MacOS.
}
  \item{infile}{
   Location of the input datafile
}
  \item{outpath}{
    location of folder to write the output files
}
  \item{numinds}{
    Number of individuals in data file
}
  \item{numloci}{
    Number of loci in data file
}
  \item{plot_output}{
 If "1" each STRUCTURE job will generate a pdf format graph in "outpath" (printqhat must be 1). If "0" no graphic output is produced
}
  \item{label}{
    Input file contains labels (names) for each individual. 1 = Yes; 0 = No
}
  \item{popdata}{
    Input file contains a user-defined population-of-origin for each individual. 1 = Yes; 0 = No
}
  \item{popflag}{
    Input file contains an indicator variable which says whether to use popinfo when USEPOPINFO==1. 1 = Yes; 0 = No
}
  \item{phenotypes}{
    Input file contains a column of phenotype information. 1 = Yes; 0 = No
}
  \item{markernames}{
    Input file contains a row of marker names. 1 = Yes; 0 = No.
}
  \item{mapdist}{
    The next row of the data file (or the first row if markernames==0) contains a list of mapdistances between neighboring loci
}
  \item{onerowperind}{
    The data for each individual are arranged in a single row. E.g., for diploid data, this would mean that the two alleles for each locus are in consecutive order in the same row, rather than being arranged in the same column, in two consecutive rows
}
  \item{phaseinfo}{
    The row(s) of genotype data for each individual are followed by a row of information about haplotype phase.
}
  \item{extracol}{
    Input file contains an extra column of data 1 = Yes; 0 = No.
}
  \item{missing}{
    code for missing genotype (Default=-9)
}
  \item{ploidy}{
    Ploidy of the organism. Default is 2 (diploid).
}
  \item{usepopinfo}{
    Use prior population information to assist clustering  1 = Yes; 0 = No.
}
  \item{printqhat}{
     the point estimate for  is not only printed into the main results file, but also into a separate file with suffix ``q''. 1 = Yes; 0 = No. This file is required in order to run the companion program STRAT and to automatically generate graphic output (see plot_output).
}

  \item{...}{
     see http://pritch.bsd.uchicago.edu/software/readme_2_1/node33.html
for complete list of arguments when running STRUCTURE from command line.
}
}
\details{
The function uses Rmpi  package to pilot serial STRUCTURE jobs in parallel  
MPI_structure requires:

1- the R package Rmpi : http://www.stats.uwo.ca/faculty/yu/Rmpi/

	- MacOS: functional version of Rmpi is available from the package repository CRAN
		(otherwise see: http://www.stats.uwo.ca/faculty/yu/Rmpi/)
		
	- Windows and LINUX see specific instruction on http://www.stats.uwo.ca/faculty/yu/Rmpi/

2- the software STRUCTURE installed: http://pritch.bsd.uchicago.edu/structure.html

}
\references{

Besnier F, Glover KA ParallelStructure: a R package to distribute parallel runs of the population genetics program STRUCTURE on multi-core computers. submitted


Rmpi website: http://www.stats.uwo.ca/faculty/yu/Rmpi/
 
 STRUCTURE website:  http://pritch.bsd.uchicago.edu/software.html  

}
\author{
Francois Besnier, Kevin Glover
}
\note{

see example data file and joblist file:

\link{structure_data}

\link{structure_jobs}


}

\examples{

## Run according to your platform (windows/Unix) and version of STRUCTURE


##FOR UNIX:

#	system('mkdir structure_results') # creat a directory to store results
#	data(structure_data)			  # call data file 
#	data(structure_jobs)			  # call joblist file

## write input files in current working directory ####

#	write(t(structure_jobs),ncol=length(structure_jobs[1,]),file='joblist1.txt')
#	write(t(structure_data),ncol=length(structure_data[1,]),file='example_data.txt')
#########################################


### call MPI_structure for the given example
### output files are stored in "structure_results/"

## You have added the path of STRUCTURE executable into /usr/local/bin:

#	MPI_structure(structure_path=NULL,joblist='joblist1.txt',n_cpu=4,infile='example_data.txt',outpath='structure_results/',numinds=987,numloci=9,printqhat=1)

##OTHERWISE

#	Mac_path="/Applications/Structure.app/Contents/Resources/Java/bin/"

## (Check the "Mac_path" variable is correct ):
#	MPI_structure(structure_path=Mac_path,joblist='joblist1.txt',n_cpu=4,infile='example_data.txt',outpath='structure_results/',numinds=987,numloci=9,printqhat=1)

##FOR WINDOWS:


#	shell('mkdir structure_results') # creat a directory to store results
#	data(structure_data)			  # call data file 
#	data(structure_jobs)			  # call joblist file

## write input files in current working directory ####
#	write(t(structure_jobs),ncol=length(structure_jobs[1,]),file='joblist1.txt')
#	write(t(structure_data),ncol=length(structure_data[1,]),file='example_data.txt')
#########################################


### call MPI_structure for the given example
### output files are stored in "structure_results/"

## You have added the path of STRUCTURE executable into your environment variable:

#	MPI_structure(structure_path=NULL,joblist='joblist1.txt',n_cpu=4,infile='example_data.txt',outpath='structure_results/',numinds=987,numloci=9,printqhat=1)

##OTHERWISE

#	Windows_path="c:/Program Files (x86)/Structure2.3.4/bin/"

##  (Check the "Windows_path" variable is correct  for your version of STRUCTURE):
#	MPI_structure(structure_path=Windows_path,joblist='joblist1.txt',n_cpu=4,infile='example_data.txt',outpath='structure_results/',numinds=987,numloci=9,printqhat=1)


}
