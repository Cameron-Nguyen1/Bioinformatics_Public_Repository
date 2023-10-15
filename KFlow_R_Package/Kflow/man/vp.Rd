\name{vp}
\alias{vp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a volcano plot using previous results from visualize and read_data
}
\description{
Makes a volcano plot using previous results alongside a comparison expression you have to manually declare. See details for more. Call ggplotly to visualize the result of this function: e.g. ggplotly(result4[[2]])
}
\usage{
vp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_from_read}{ = Reference to object obtained from "$read_data()" call

}
  \item{data_from_vis}{ = Reference to object obtained from "$visualize()" call

}
  \item{labels}{ = Item [[3]] from "$visualize()" call, e.g. visualize_object[[3]]

}
}
\details{ = This function requires you to declare a variable named "ex" before calling the function.
This variable is an expression used to compare one category against the other, typically a comparison of test vs control. "infection" is mandatory in this expression.
e.g. ex = "infection=HI7_100nM_72hpi-IFNAR1_100nM_72hpi"
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
\author{
Cameron Nguyen
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
result3 = k$pca(result1, result2) #Uses result2 return value, influenced by option 1, 2, or 3 picked in k$visualize() call
result3 #Call result object to get map
ex = "infection=HI7_100nM_72hpi-IFNAR1_100nM_72hpi" #Define expression for Volcano Plot **Group** Comparison // Not sample, group!
result4 = k$vp(result1, result2, result2[[3]]) # volcano plot using result from visualize method, uses ex to define comparison
ggplotly(result4[[2]]) #Plot the Volcano plot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
