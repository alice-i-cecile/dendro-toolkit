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
standardize_mle(tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", error = "lnorm", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
}
  \item{error}{
%%     ~~Describe \code{error} here~~
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
function (tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", 
    error = "lnorm", ...) 
{
    dummy_effects <- pad_effects(list(), tra, form)
    flat_effects0 <- flatten_effects(dummy_effects, model)
    fes_likelihood <- function(flat_effects) {
        effects <- unflatten_effects(flat_effects)
        effects <- pad_effects(effects, tra, form)
        predicted <- predicted_tra(effects, tra, form)
        residuals <- residuals_tra(tra, predicted, error)
        llh <- llh_tra(residuals, error)
        return(-llh)
    }
    mle_solution <- optim(flat_effects0, fes_likelihood, ...)
    print(mle_solution[2:5])
    effects <- unflatten_effects(mle_solution$par)
    effects <- pad_effects(effects, tra, form)
    effects <- sort_effects(effects, tra)
    effects <- rescale_effects(effects, form)
    fit <- model_fit_tra(effects, tra, model, form, error)
    settings <- list(model = model, form = form, error = error, 
        method = "likelihood", ...)
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
