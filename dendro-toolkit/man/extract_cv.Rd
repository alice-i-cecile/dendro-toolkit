\name{extract_cv}
\alias{extract_cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_cv(growth_model, model_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growth_model}{
%%     ~~Describe \code{growth_model} here~~
}
  \item{model_type}{
%%     ~~Describe \code{model_type} here~~
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
function (growth_model, model_type) 
{
    model_class <- class(growth_model)
    if (model_type == "gam") {
        cv <- extract_cv_gam(growth_model)
    }
    else if (model_type == "gnm") {
        cv <- extract_cv_gnm(growth_model)
    }
    else if (model_type == "lm") {
        cv <- extract_cv_lm(growth_model)
    }
    return(cv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
