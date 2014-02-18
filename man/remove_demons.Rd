\name{remove_demons}
\alias{remove_demons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
remove_demons(m = 1, corrupted_effects, tra, sparse = TRUE, form = "multiplicative")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{corrupted_effects}{
%%     ~~Describe \code{corrupted_effects} here~~
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
function (m = 1, corrupted_effects, tra, sparse = TRUE, form = "multiplicative") 
{
    ridge_ids <- demonic_ritual(corrupted_effects, tra, sparse)
    effects <- corrupted_effects
    if (form == "multiplicative") {
        effects <- lapply(effects, log)
    }
    effects$I <- effects$I - m * ridge_ids$I
    effects$T <- effects$T - m * ridge_ids$T
    effects$A <- effects$A + m * ridge_ids$A
    if (form == "multiplicative") {
        effects <- lapply(effects, exp)
    }
    effects <- rescale_effects(effects, form)
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
