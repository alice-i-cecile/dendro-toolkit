\name{extract_cv_gam}
\alias{extract_cv_gam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_cv_gam(growth_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growth_model}{
%%     ~~Describe \code{growth_model} here~~
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
function (growth_model) 
{
    coeff <- dummy.coef(growth_model)
    Q.cv <- coeff$Q.cv
    F.cv <- coeff$F.cv
    tab.tra <- growth_model$data
    smooth_fit <- growth_model$smooth
    unique_levels <- unique(tab.tra$A.num)
    indices <- unlist(foreach(i = unique_levels) \%do\% {
        match(i, tab.tra$A.num)
    })
    A.cv <- smooth_fit[indices]
    names(A.cv) <- unique_levels
    cv <- list(Q = Q.cv, F = F.cv, A = A.cv)
    return(cv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
