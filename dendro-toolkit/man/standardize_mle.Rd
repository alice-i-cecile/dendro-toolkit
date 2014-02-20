\name{standardize_mle}
\alias{standardize_mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_mle(tra, model = c("Time", "Age"), link = "log", method = "CG", ...)
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
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (tra, model = c("Time", "Age"), link = "log", method = "CG", 
    ...) 
{
    effects0 <- vector(mode = "list", length = length(model))
    names(effects) <- model
    for (i in model) {
        dim_i <- nlevels(tra[[i]])
        if (form == "additive") {
            effects[[i]] <- rep.int(0, dim_i)
        }
        else {
            effects[[i]] <- rep.int(1, dim_i)
        }
        names(effects[[i]]) <- levels(tra[[i]])
    }
    fes_likelihood <- function(effects) {
        predicted <- predicted_tra(effects, tra, link)
        residuals <- residuals_tra(tra, predicted, link)
        llh <- llh_tra(residuals, error)
        return(llh)
    }
    mle_solution <- optim(effects0, -fes_likelihood, ...)
    print(mle_solution[2:5])
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
