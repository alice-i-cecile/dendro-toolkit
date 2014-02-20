\name{estimate_demons}
\alias{estimate_demons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_demons(corrupted_effects, tra = TRUE, link = "log")
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
    ridge_ids <- demonic_ritual(corrupted_effects, tra)
    if (link == "log") {
        corrupted_effects <- lapply(corrupted_effects, log)
    }
    demonic_llh <- function(purified_effects) {
        sd_I <- sd(purified_effects$Timeree)
        sd_T <- sd(purified_effects$Time)
        sd_A <- sd(purified_effects$Age)
        llh_I <- sum(dnorm(purified_effects$Timeree, sd = sd_I, 
            log = TRUE))
        llh_T <- sum(dnorm(purified_effects$Time, sd = sd_T, 
            log = TRUE))
        llh_A <- sum(dnorm(purified_effects$Age, mean = mean(purified_effects$Age), 
            sd = sd_A, log = TRUE))
        llh <- sum(llh_I, llh_T, llh_A)
        return(llh)
    }
    purify_probe <- function(m) {
        purified_effects <- remove_demons(m, corrupted_effects, 
            tra, link = "additive")
        llh <- demonic_llh(purified_effects)
        return(-llh)
    }
    optimal_m <- optimize(f = purify_probe, interval = c(-1, 
        1))$minimum
    return(optimal_m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
