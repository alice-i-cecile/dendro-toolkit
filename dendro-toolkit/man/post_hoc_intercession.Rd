\name{post_hoc_intercession}
\alias{post_hoc_intercession}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
post_hoc_intercession(corrupted_effects, tra = TRUE, link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corrupted_effects}{
%%     ~~Describe \code{corrupted_effects} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
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
function (corrupted_effects, tra = TRUE, link = "log") 
{
    optimal_m <- estimate_demons(corrupted_effects, tra, link)
    if (abs(optimal_m) > 0.99) {
        return(corrupted_effects)
    }
    print(paste("m of", optimal_m, "selected in choosing solution on likelihood ridge."))
    purified <- remove_demons(optimal_m, corrupted_effects, tra, 
        link)
    return(purified)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
