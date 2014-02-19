\name{standardize_tsfs}
\alias{standardize_tsfs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_tsfs(tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", error = "lnorm", cor_threshold = 0.999999)
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
  \item{cor_threshold}{
%%     ~~Describe \code{cor_threshold} here~~
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
    error = "lnorm", cor_threshold = 0.999999) 
{
    if (error == "lnorm") {
        mean_type <- "geometric"
    }
    else {
        mean_type <- "arithmetic"
    }
    effects <- vector(mode = "list", length = 3)
    names(effects) <- c("I", "T", "A")
    for (i in 1:3) {
        dim_i <- nlevels(tra[[i + 1]])
        if (form == "additive") {
            effects[[i]] <- rep.int(0, dim_i)
        }
        else {
            effects[[i]] <- rep.int(1, dim_i)
        }
        names(effects[[i]]) <- levels(tra[[i + 1]])
    }
    inc_effects <- names(model[unlist(model) == TRUE])
    name_dim_dict <- c(I = "I", T = "T", A = "A")
    effect_order <- match(inc_effects, name_dim_dict)
    if (length(inc_effects) == 0) {
        converged <- TRUE
    }
    else {
        converged <- FALSE
    }
    iteration <- 0
    working_tra <- tra
    while (!converged) {
        last_tra <- working_tra
        iteration <- iteration + 1
        print(paste("Iteration", iteration))
        for (j in effect_order) {
            est_j <- est_effect(working_tra, j, mean_type)
            working_tra <- remove_effect(working_tra, est_j, 
                j, form)
            if (form == "additive") {
                effects[[j]] <- effects[[j]] + est_j
            }
            else {
                effects[[j]] <- effects[[j]] * est_j
            }
        }
        if (error == "norm") {
            conv_cor <- cor(working_tra$G, last_tra$G)
        }
        else {
            conv_cor <- cor(log(working_tra$G), log(last_tra$G))
        }
        if (conv_cor >= cor_threshold) {
            converged <- TRUE
        }
    }
    effects <- sort_effects(effects, tra)
    effects <- rescale_effects(effects, form)
    fit <- model_fit_tra(effects, tra, model, form, error)
    settings <- list(model = model, form = form, error = error, 
        method = "sfs")
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
