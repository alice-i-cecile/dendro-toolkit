\name{flatten_effects}
\alias{flatten_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
flatten_effects(effects, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (effects, model) 
{
    flat_effects <- vector()
    if (model$I) {
        effects_I <- dummy_effects$I
        names(effects_I) <- paste("I", names(dummy_effects$I), 
            sep = "@")
        flat_effects <- c(flat_effects, effects_I)
    }
    if (model$T) {
        effects_T <- dummy_effects$T
        names(effects_T) <- paste("T", names(dummy_effects$T), 
            sep = "@")
        flat_effects <- c(flat_effects, effects_T)
    }
    if (model$A) {
        effects_A <- dummy_effects$A
        names(effects_A) <- paste("A", names(dummy_effects$A), 
            sep = "@")
        flat_effects <- c(flat_effects, effects_A)
    }
    return(flat_effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
