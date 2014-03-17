\name{grab_age}
\alias{grab_age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the age of the tree at time of sampling
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grab_age(i, tra, sparse = TRUE)
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
        age <- max(as.numeric(as.character(tra[tra$i == i, "a"])))
    }
    else {
        age_row <- min(which(!is.na(tra[i, , ]), arr.ind = T)[[2]])
        age <- dimnames(tra)[[3]][age_row]
    }
    return(as.numeric(age))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
