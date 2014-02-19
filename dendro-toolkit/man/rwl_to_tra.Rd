\name{rwl_to_tra}
\alias{rwl_to_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rwl_to_tra(rwl, birth_years = NULL)
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
    filled <- which(!is.na(rwl), arr.ind = TRUE)
    n_data_points <- nrow(filled)
    add_growth_data <- function(year, tree) {
        birth <- birth_index[tree]
        G <- rwl[year, tree]
        i <- names(rwl)[tree]
        t <- rownames(rwl)[year]
        a <- year - birth + 1
        out <- list(G, i, t, a)
        return(out)
    }
    raw_tra <- mapply(add_growth_data, year = filled[, 1], tree = filled[, 
        2], SIMPLIFY = FALSE)
    tra <- data.frame(matrix(unlist(raw_tra), ncol = 4, byrow = TRUE))
    names(tra) <- c("Growth", "Tree", "Year", "Age")
    tra$Growth <- as.numeric(tra$Growth)
    return(tra)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
