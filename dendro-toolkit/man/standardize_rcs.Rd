\name{standardize_rcs}
\alias{standardize_rcs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_rcs(tra, model = list(I = FALSE, A = TRUE, T = TRUE), form = "multiplicative", error = "lnorm", sparse = TRUE)
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
function (tra, model = list(I = FALSE, A = TRUE, T = TRUE), form = "multiplicative", 
    error = "lnorm", sparse = TRUE) 
{
    if (sparse) {
        if (!is.data.frame(tra)) {
            tra <- sparse_tra(tra)
        }
    }
    else {
        if (is.data.frame(tra)) {
            tra <- unsparse_tra(tra)
        }
    }
    if (error == "lnorm") {
        mean_type <- "geometric"
    }
    else {
        mean_type <- "arithmetic"
    }
    inc_effects <- names(model[unlist(model) == TRUE])
    name_dim_dict <- c(I = "I", T = "T", A = "A")
    effect_order <- match(inc_effects, name_dim_dict)
    effects <- list(I = NULL, T = NULL, A = NULL)
    working_tra <- tra
    for (effect in effect_order) {
        effects[[effect]] <- est_effect(working_tra, effect, 
            mean_type, sparse)
        working_tra <- remove_effect(working_tra, effects[[effect]], 
            effect, form, sparse)
    }
    effects <- pad_effects(effects, tra, form, sparse)
    effects <- sort_effects(effects, tra, sparse)
    effects <- rescale_effects(effects, form)
    fit <- model_fit_tra(effects, tra, model, form, error, sparse)
    settings <- list(model = model, form = form, error = error, 
        sparse = sparse, method = "rcs")
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
