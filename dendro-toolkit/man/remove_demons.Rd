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
remove_demons(m = 1, corrupted_effects, tra = TRUE, link = "log")
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
function (m = 1, corrupted_effects, tra = TRUE, link = "log") 
{
    ridge_ids <- demonic_ritual(corrupted_effects, tra)
    effects <- corrupted_effects
    if (link == "log") {
        effects <- lapply(effects, log)
    }
    effects$Timeree <- effects$Timeree - m * ridge_ids$Timeree
    effects$Time <- effects$Time - m * ridge_ids$Time
    effects$Age <- effects$Age + m * ridge_ids$Age
    if (link == "log") {
        effects <- lapply(effects, exp)
    }
    effects <- rescale_effects(effects, link)
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
