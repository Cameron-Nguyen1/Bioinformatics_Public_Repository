\name{plot_switch}
\alias{plot_switch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots your switch graphs!
}
\description{ plot_switch accepts the return object of DTU_analysis which is a switch list. The output includes switch graphs with predicted consequences. }
\usage{
plot_switch(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{switch_list}{ = A switch list object returned from DTU_Analysis() call.
}
  \item{gene_name}{ = The gene name you want the switch graph for. See "Switch_top_50.csv" for genes.
}
  \item{plot_name}{ = The name of the image file you'll save after calling plot_switch().
}
  \item{cond1}{ = condition 1, i.e. healthy / unhealthy
}
  \item{cond2}{ = condition 2.
}
}
\details{ Plot the DTU switch plots!
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
k$plot_switch(myswitch, 'GABRD', "GABRDSwitchPlot.png", 'Treated','Untreated')
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
