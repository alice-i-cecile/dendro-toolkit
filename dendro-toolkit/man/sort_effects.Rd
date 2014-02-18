\name{sort_effects}
\alias{sort_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sort_effects(effects, tra, sparse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
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
function (effects, tra, sparse = TRUE) 
{
    if (sparse) {
        sorted_effects <- list()
        effect_names <- names(effects$I)
        sorted_effects$I <- effects$I[sort(effect_names)]
        for (j in c("T", "A")) {
            effect_names <- as.numeric(names(effects[[j]]))
            sorted_effects[[j]] <- effects[[j]][as.character(sort(effect_names))]
        }
    }
    else {
        sorted_effects <- effects
        for (i in 1:3) {
            correct_order <- dimnames(tra)[[i]]
            sorted_effects[[i]] <- effects[[i]][correct_order]
        }
    }
    return(sorted_effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
