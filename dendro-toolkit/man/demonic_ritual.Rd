\name{demonic_ritual}
\alias{demonic_ritual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
demonic_ritual(effects, tra = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
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
function (effects, tra = TRUE) 
{
    birth_years <- get_birth_years(tra)
    birth_index <- sapply(birth_years, get_birth_index, tra = tra)
    k <- -birth_index
    years <- sort(as.numeric(levels(tra$t)))
    year_index <- match(names(effects$T), as.character(years))
    names(year_index) <- names(effects$T)
    ages <- sort(as.numeric(levels(tra$a)))
    age_index <- match(names(effects$A), as.character(ages))
    names(age_index) <- names(effects$A)
    k <- k[!is.na(k)]
    year_index <- year_index[!is.na(year_index)]
    age_index <- age_index[!is.na(age_index)]
    ridge_ids <- list(I = k, T = year_index, A = age_index)
    return(ridge_ids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
