\name{remove_effect}
\alias{remove_effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
remove_effect(tra, effect, factor.dim, form = "multiplicative", sparse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{effect}{
%%     ~~Describe \code{effect} here~~
}
  \item{factor.dim}{
%%     ~~Describe \code{factor.dim} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
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
function (tra, effect, factor.dim, form = "multiplicative", sparse = TRUE) 
{
    if (sparse) {
        id <- factor.dim + 1
        removed_tra <- tra
        for (effect_id in names(effect)) {
            relevant_rows <- tra[[id]] == effect_id
            if (form == "additive") {
                removed_tra[relevant_rows, "G"] <- removed_tra[relevant_rows, 
                  "G"] - effect[effect_id]
            }
            else {
                removed_tra[relevant_rows, "G"] <- removed_tra[relevant_rows, 
                  "G"]/effect[effect_id]
            }
        }
    }
    else {
        FUN <- ifelse(form == "additive", "-", "/")
        removed_tra <- sweep(tra, factor.dim, effect, FUN)
    }
    return(removed_tra)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
