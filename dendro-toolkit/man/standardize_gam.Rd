\name{standardize_gam}
\alias{standardize_gam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_gam(tra, model = c("Time", "Age"), link = "log", sparse = TRUE, age_k = 10, ...)
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
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
}
  \item{age_k}{
%%     ~~Describe \code{age_k} here~~
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
function (tra, model = c("Time", "Age"), link = "log", sparse = TRUE, 
    age_k = 10, ...) 
{
    tra$Age <- as.numeric(as.character(tra$Age))
    growth_formula <- as.formula(make_gam_formula(model))
    print("Model initialized.")
    print(growth_formula)
    if (link == "identity") {
        family <- gaussian(link = "identity")
    }
    else if (link == "log") {
        tra$Growth <- log(tra$Growth)
        family <- gaussian(link = "identity")
    }
    growth_model <- gam(growth_formula, family = family, data = tra, 
        ...)
    print("Generalized additive model used to standardize data.")
    effects <- extract_effects_gam(growth_model, model, link == 
        "log", tra)
    if (form == "multiplicative") {
        effects <- lapply(effects, exp)
        tra$Growth <- exp(tra$Growth)
    }
    print("Effects extracted.")
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
