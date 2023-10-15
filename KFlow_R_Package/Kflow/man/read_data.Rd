\name{read_data}
\alias{read_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in Kallisto data
%%  ~~function to do ... ~~
}
\description{ Reads in Kallisto data, enables Kflow to execute next steps
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_data(wdir, study_file, db)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wdir}{
     = Working Directory with Kallisto Quant folders
}
  \item{study_file}{
     = Study Design File
}
  \item{db}{
     = Database Object, tested with Ensembl database from R library, see Examples below
  }
}
\details{
%%
}
\value{
%% NULL
}
\references{
%% Berry, Alexander S., et al. “An Open-Source Toolkit to Expand Bioinformatics Training in Infectious Diseases.” MBio, vol. 12, no. 4, 2021, https://doi.org/10.1128/mbio.01214-21.
}
\author{
%%  Cameron Nguyen
}
\note{
%%
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~
}
\examples{
#result1 = k$read_data("/Users/azn/Documents/HAE_quants_compare/HI7_100_72", "Study_Design_Redo.tsv", EnsDb.Hsapiens.v86)

## The function is currently defined as
function (x)
{
  }
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
