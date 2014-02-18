\name{sample_depth_tra}
\alias{sample_depth_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample_depth_tra(tra, factor.dim = 2, sparse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{factor.dim}{
%%     ~~Describe \code{factor.dim} here~~
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
function (tra, factor.dim = 2, sparse = FALSE) 
{
    if (sparse) {
        if (!is.data.frame(tra)) {
            tra <- sparse_tra(tra)
        }
    }
    else {
        if (is.data.frame(tra)) {
            tra <- unsparse_tra(tra)
        }
    }
    if (sparse) {
        positions <- tra[[factor.dim + 1]]
        ids <- levels(positions)
        sample_depth <- sapply(ids, function(x) {
            sum(positions == x)
        })
        if (factor.dim == 1) {
            ordering <- sort(names(sample_depth))
        }
        else {
            ordering <- sort(as.numeric(names(sample_depth)))
        }
        sample_depth <- sample_depth[sapply(ordering, as.character)]
    }
    else {
        filled_cells <- !is.na(tra)
        sample_depth <- apply(filled_cells, factor.dim, sum)
    }
    return(sample_depth)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
