\name{grab_birth_year}
\alias{grab_birth_year}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grab_birth_year(i, tra, sparse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
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
function (i, tra, sparse = TRUE) 
{
    if (sparse) {
        slice <- tra[tra$i == i, ][1, ]
        birth_year <- as.numeric(as.character(slice$t)) - as.numeric(as.character(slice$a)) + 
            1
    }
    else {
        slice <- which(!is.na(tra[i, , ]), arr.ind = T)[1, ]
        birth_index <- slice[1] - slice[2] + 1
        birth_year <- as.numeric(dimnames(tra)[[2]][birth_index])
    }
    return(birth_year)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
