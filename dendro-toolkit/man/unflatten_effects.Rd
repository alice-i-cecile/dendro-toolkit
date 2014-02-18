\name{unflatten_effects}
\alias{unflatten_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unflatten_effects(flat_effects)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flat_effects}{
%%     ~~Describe \code{flat_effects} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (flat_effects) 
{
    effects <- list()
    id_code <- substr(names(flat_effects), 1, 1)
    names(flat_effects) <- substr(names(flat_effects), 3, length(names(flat_effects)))
    effects$I <- flat_effects[id_code == "I"]
    effects$T <- flat_effects[id_code == "T"]
    effects$A <- flat_effects[id_code == "A"]
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
