\name{plot_str}
\alias{plot_str}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PLOT FUNCTION FOR STRUCTURE OUTPUT
}
\description{
internal function of MPI_structure. if "plot_output" argument of MOI_structure is 1, "plot_stru()" is called to generate structure plots.
}
\usage{
plot_str(loc = NULL, list_jobs, x_width, y_height)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
Location to write graphic output
}
  \item{list_jobs}{
list of jobs from which output will be plotted
}
  \item{x_width}{
width of the graph
}
  \item{y_height}{
height of the graph
}
}
\author{
Francois Besnier, Kevin Glover
}
