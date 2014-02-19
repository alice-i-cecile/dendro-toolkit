\name{tra_to_rwl}
\alias{tra_to_rwl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tra_to_rwl(tra, data_column = "Growth")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{data_column}{
%%     ~~Describe \code{data_column} here~~
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
function (tra, data_column = "Growth") 
{
    tree_names <- unique(tra$Tree)
    years <- unique(tra$Year)
    years <- sort(as.numeric(as.character(year)))
    num_trees <- length(tree_names)
    num_years <- length(years)
    rwl <- as.data.frame(matrix(data = NA, nrow = num_year, ncol = num_trees, 
        dimnames = list(years, trees)))
    for (i in 1:nrows(tra)) {
        tree <- tra[i, "Tree"]
        year <- tra[i, "Year"]
        tree_index <- which(names(rwl) == tree)
        year_index <- which(rownames(rwl) == year)
        rwl[year_index, tree_index] <- tra[i, data_column]
    }
    return(rwl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
