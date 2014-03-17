\name{tra_to_rwl}
\alias{tra_to_rwl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts .tra format tree ring data to .rwl format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tra.to.rwl(tra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
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
function (tra) 
{
    i.size <- dim(tra)[1]
    t.size <- dim(tra)[2]
    rwl <- as.data.frame(matrix(NA, t.size, i.size))
    colnames(rwl) <- dimnames(tra)[[1]]
    rownames(rwl) <- dimnames(tra)[[2]]
    for (tree in 1:i.size) {
        filled <- !is.na(tra[tree, , ])
        filledYears <- rownames(filled[rowSums(filled) > 0, ])
        data <- tra[tree, , ]
        data <- data[!is.na(data)]
        rwl[filledYears, tree] <- data
    }
    return(rwl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
