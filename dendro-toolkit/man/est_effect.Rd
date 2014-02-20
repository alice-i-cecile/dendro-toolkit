\name{est_effect}
\alias{est_effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
est_effect(tra, id, mean_type = "arithmetic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{mean_type}{
%%     ~~Describe \code{mean_type} here~~
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
function (tra, id, mean_type = "arithmetic") 
{
    estimate_effect <- function(id_i) {
        data <- tra[tra[[id]] == id_i, "G"]
        if (mean_type == "geometric") {
            est_effect <- geomMean(data)
        }
        else {
            est_effect <- mean(data, na.rm = TRUE)
        }
        return(est_effect)
    }
    est_effect <- sapply(levels(tra[[id]]), estimate_effect)
    return(est_effect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
