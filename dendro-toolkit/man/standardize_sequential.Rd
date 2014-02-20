\name{standardize_sequential}
\alias{standardize_sequential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_sequential(tra, model = c("Age", "Time"), link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
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
function (tra, model = c("Age", "Time"), link = "log") 
{
    if (link == "log") {
        error <- "lnorm"
        form <- "multiplicative"
    }
    else if (link == "identity") {
        error <- "norm"
        form <- "additive"
    }
    if (error == "lnorm") {
        mean_type <- "geometric"
    }
    else {
        mean_type <- "arithmetic"
    }
    inc_effects <- model
    effects <- vector(mode = "list", length = length(model))
    names(effects) <- model
    working_tra <- tra
    for (id in model) {
        effects[[id]] <- est_effect(working_tra, id, mean_type)
        working_tra <- remove_effect(working_tra, effects[[id]], 
            effect, link)
    }
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
