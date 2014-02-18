\name{unsparse_tra}
\alias{unsparse_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unsparse_tra(stra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stra}{
%%     ~~Describe \code{stra} here~~
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
function (stra) 
{
    indices <- lapply(stra[2:ncol(stra)], levels)
    indices[c("t", "a")] <- lapply(indices[c("t", "a")], as.numeric)
    indices <- lapply(indices, sort)
    indices[c("t", "a")] <- lapply(indices[c("t", "a")], as.character)
    tra <- array(data = NA, dim = sapply(indices, length), dimnames = indices)
    for (j in 1:nrow(stra)) {
        stra_row <- stra[j, ]
        i <- as.character(stra[j, "i"])
        t <- as.character(stra[j, "t"])
        a <- as.character(stra[j, "a"])
        G <- stra[j, "G"]
        tra[i, t, a] <- G
    }
    return(tra)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
