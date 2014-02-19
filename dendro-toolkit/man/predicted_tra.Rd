\name{predicted_tra}
\alias{predicted_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predicted_tra(effects, tra, form)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
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
function (effects, tra, form) 
{
    predicted <- tra
    for (r in 1:nrow(predicted)) {
        i <- as.character(tra[r, "i"])
        t <- as.character(tra[r, "t"])
        a <- as.character(tra[r, "a"])
        if (form == "additive") {
            predicted[r, "G"] <- effects$I[i] + effects$T[t] + 
                effects$A[a]
        }
        else {
            predicted[r, "G"] <- effects$I[i] * effects$T[t] * 
                effects$A[a]
        }
    }
    return(predicted)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
