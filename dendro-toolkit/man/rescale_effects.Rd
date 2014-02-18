\name{rescale_effects}
\alias{rescale_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rescale_effects(effects, form = "multiplicative")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
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
function (effects, form = "multiplicative") 
{
    if (form == "multiplicative") {
        effects <- lapply(effects, log)
    }
    mean_effects <- lapply(effects, mean, na.rm = T)
    if (!sum(!is.null(effects[[3]]))) {
        return(effects)
    }
    if (sum(!is.null(effects[[1]]))) {
        effects[[1]] <- effects[[1]] - mean_effects[[1]]
        effects[[3]] <- effects[[3]] + mean_effects[[1]]
    }
    if (sum(!is.null(effects[[2]]))) {
        effects[[2]] <- effects[[2]] - mean_effects[[2]]
        effects[[3]] <- effects[[3]] + mean_effects[[2]]
    }
    if (form == "multiplicative") {
        effects <- lapply(effects, exp)
    }
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
