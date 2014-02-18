\name{add_demons}
\alias{add_demons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add_demons(m = 1, effects, tra, sparse = TRUE, form = "multiplicative")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
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
function (m = 1, effects, tra, sparse = TRUE, form = "multiplicative") 
{
    ridge_ids <- demonic_ritual(effects, tra, sparse)
    corrupted_effects <- effects
    if (form == "multiplicative") {
        corrupted_effects <- lapply(corrupted_effects, log)
    }
    corrupted_effects$I <- corrupted_effects$I + m * ridge_ids$I
    corrupted_effects$T <- corrupted_effects$T + m * ridge_ids$T
    corrupted_effects$A <- corrupted_effects$A - m * ridge_ids$A
    if (form == "multiplicative") {
        corrupted_effects <- lapply(corrupted_effects, exp)
    }
    corrupted_effects <- rescale_effects(corrupted_effects, form)
    return(corrupted_effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
