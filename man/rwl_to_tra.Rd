\name{rwl_to_tra}
\alias{rwl_to_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts .rwl format tree ring data to .tra format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rwl.to.tra(rwl, birth_years = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rwl}{
%%     ~~Describe \code{rwl} here~~
}
  \item{birth_years}{
%%     ~~Describe \code{birth_years} here~~
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
function (rwl, birth_years = NULL) 
{
    if (is.null(birth_years)) {
        birth_years <- foreach(i = colnames(rwl)) \%do\% {
            getEndpoints(rwl[i], side = "start")
        }
        names(birth_years) <- names(rwl)
    }
    else {
        birth_limit <- min(sapply(birth_years, as.numeric))
        rwl_limit <- min(sapply(rownames(rwl), as.numeric))
        if (birth_limit < rwl_limit) {
            num_missing <- rwl_limit - birth_limit
            empty_rows <- rbind(rwl[0, ], matrix(NA, num_missing, 
                ncol(rwl)))
            colnames(empty_rows) <- colnames(rwl)
            rownames(empty_rows) <- birth_limit:(rwl_limit - 
                1)
            rwl <- rbind(empty_rows, rwl)
        }
    }
    birth_index <- vector()
    for (i in 1:length(birth_years)) {
        birth_index[i] <- which(rownames(rwl) == toString(birth_years[i]))
    }
    names(birth_index) <- names(birth_years)
    i.size <- ncol(rwl)
    t.size <- nrow(rwl)
    T.size <- t.size
    tra <- array(NA, c(i.size, t.size, T.size), list(colnames(rwl), 
        rownames(rwl), 1:T.size))
    for (tree in 1:ncol(rwl)) {
        birth <- birth_index[tree]
        for (year in 1:nrow(rwl)) {
            age <- year - birth + 1
            datum <- rwl[year, tree]
            if (!is.na(datum)) {
                tra[tree, year, age] <- datum
            }
        }
    }
    empty_tra <- is.na(tra)
    empty_Q <- apply(empty_tra, 1, all)
    empty_F <- apply(empty_tra, 2, all)
    empty_A <- apply(empty_tra, 3, all)
    tra <- tra[!empty_Q, !empty_F, !empty_A]
    return(tra)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
