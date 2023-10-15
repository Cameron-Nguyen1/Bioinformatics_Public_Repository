\name{visualize}
\alias{visualize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize data after reading it in
}
\description{
visualize aims to normalize and visualize quantification data. Must be used after "read_data".
}
\usage{
visualize(read_res, col_start, col_end, graph_option)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read_res}{ = Result object from "$read_data()" call.

}
  \item{col_start}{ = Beginning column corresponding to the first sample in your Study Design File.

}
  \item{col_end}{ = Ending column corresponding to the final sample in your Study Design File.

}
  \item{graph_option}{ = Options of int 1, 2, or 3. 1 = Unfiltered, nonnormalized graph of CPM
  2 = Filtered nonnormalized graph of CPM
  3 = Filited, normalzied graph of CPM

}
}
\details{
%% Visualize your data from the step 1 "read_data" call.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% Berry, Alexander S., et al. “An Open-Source Toolkit to Expand Bioinformatics Training in Infectious Diseases.” MBio, vol. 12, no. 4, 2021, https://doi.org/10.1128/mbio.01214-21.
}
\author{Cameron Nguyen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
k = kalli_flow() #Instantiate the class
result1 = k$read_data("/Users/azn/Documents/HAE_quants_compare/HI7_100_72", "Study_Design_Redo.tsv", EnsDb.Hsapiens.v86) #Read data (Args=WorkingDir, Design_File, Database for annotation)
result2 = k$visualize(result1, "HI7_100nM_72hpi_85", "IFNAR1_100nM_72hpi_99",3) #Visualize results ARGS=read_data obj, col_start, col_end, (1, 2 or 3) 1 = Non-norm, unfilt, 2 = non-norm, filtered, 3 = norm, filtered
}
