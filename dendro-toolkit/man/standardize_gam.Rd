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
standardize_gam(tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", error = "lnorm", sparse = TRUE, age_k = 10, ...)
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
function (tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", 
    error = "lnorm", sparse = TRUE, age_k = 10, ...) 
{
    if ((form == "additive" & error == "lnorm") | (form == "multiplicative" & 
        error == "norm")) {
        simpleError("Model form and error distribution are an unrealistic match. Generalized additive models cannot be used for this configuration.")
    }
    if (sparse) {
        if (!is.data.frame(tra)) {
            tra <- sparse_tra(tra)
        }
    }
    else {
        simpleError("GAM standardization can only handle sparse tree-ring arrays.")
    }
    growth_formula <- as.formula(make_gam_formula(model))
    print("Model initialized.")
    print(growth_formula)
    if (form == "additive") {
        family <- gaussian(link = "identity")
    }
    if (form == "multiplicative") {
        tra$G <- log(tra$G)
        family <- gaussian(link = "identity")
    }
    growth_model <- gam(growth_formula, family = family, data = tra, 
        ...)
    print("Generalized additive model used to standardize data.")
    effects <- extract_effects_gam(growth_model, model, form = "additive", 
        tra)
    if (form == "multiplicative") {
        effects <- lapply(effects, exp)
        tra$G <- exp(tra$G)
    }
    print("Effects extracted.")
    effects <- pad_effects(effects, tra, form, sparse)
    effects <- sort_effects(effects, tra, sparse)
    effects <- rescale_effects(effects, form)
    k <- sum(growth_model$edf) + 1
    fit <- model_fit_tra(effects, tra, model, form, error, sparse, 
        method = "gam", k = k)
    print("Fit computed.")
    settings <- list(model = model, form = form, error = error, 
        sparse = sparse, method = "gam", ...)
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
