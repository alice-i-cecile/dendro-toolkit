\name{standardize_tra}
\alias{standardize_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize tree-ring data and produce a chronology
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_tra(tra, model = c("Age", "Time"), link = "log", optim = "sfs", post_hoc = TRUE, return_data = FALSE, make_plots = TRUE, ...)
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
  \item{optim}{
%%     ~~Describe \code{optim} here~~
}
  \item{post_hoc}{
%%     ~~Describe \code{post_hoc} here~~
}
  \item{return_data}{
%%     ~~Describe \code{return_data} here~~
}
  \item{make_plots}{
%%     ~~Describe \code{make_plots} here~~
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
function (tra, model = c("Age", "Time"), link = "log", optim = "sfs", 
    post_hoc = TRUE, return_data = FALSE, make_plots = TRUE, 
    ...) 
{
    if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra <= 
        0], na.rm = TRUE) > 0)) {
        if (link == "log") {
            stop("Zero or negative values cannot be use. Estimated effects will not be stable.")
        }
    }
    if (optim == "likelihood") {
        effects <- standardize_likelihood(tra, model, link, ...)
    }
    else if (optim == "alternate") {
        effects <- standardize_alternate(tra, model, link, ...)
    }
    else if (optim == "rcs") {
        effects <- standardize_sequential(tra, model, link, ...)
    }
    else if (optim == "gam") {
        effects <- standardize_gam(tra, model, link, ...)
    }
    if (("Tree" \%in\% model) & ("Time" \%in\% model) & ("Age" \%in\% 
        model)) {
        if (post_hoc) {
            effects <- post_hoc_intercession(out$effects, out$tra, 
                link)
            warning("Tree-time-age model selected. Post-hoc effect selection was used to stabilize parameter estimates.")
        }
        else {
            warning("Tree-time-age model selected. Parameter estimates will be unreliable. Consider using post-hoc effect selection.")
        }
    }
    effects <- sort_effects(effects, tra, sparse)
    effects <- rescale_effects(effects, link)
    fit <- model_fit_tra(effects, tra, model, link)
    data <- list(original = tra, predicted = fit$predicted, residuals = fit$residuals)
    fit <- fit[-which(names(fit) == "predicted" | names(fit) == 
        "residuals")]
    if (make_plots) {
        plots <- make_standardization_plots(effects, data, link)
        print(plots)
    }
    settings <- list(model = model, link = link, optim = optim)
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    if (return_data) {
        out <- c(out, list(data = data))
    }
    if (make_plots) {
        out <- c(out, list(plots = plots))
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
