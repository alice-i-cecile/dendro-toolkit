\name{standardize_sfs}
\alias{standardize_sfs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_sfs(tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", error = "lnorm", sparse = TRUE, cor_threshold = 0.999999)
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
    error = "lnorm", sparse = TRUE, cor_threshold = 0.999999) 
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
    converged <- FALSE
    iteration <- 0
    working_tra <- tra
    while (!converged) {
        last_tra <- working_tra
        iteration <- iteration + 1
        print(paste("Iteration", iteration))
        est_1 <- est_effect(working_tra, effect_order[1], mean_type, 
            sparse)
        intermediate_tra <- remove_effect(working_tra, est_1, 
            effect_order[1], form, sparse)
        est_2 <- est_effect(intermediate_tra, effect_order[2], 
            mean_type, sparse)
        working_tra <- remove_effect(working_tra, est_2, effect_order[2], 
            form, sparse)
        if (error == "norm") {
            if (sparse) {
                conv_cor <- cor(working_tra$G, last_tra$G)
            }
            else {
                conv_cor <- cor(working_tra, last_tra, "complete.obs")
            }
        }
        else {
            if (sparse) {
                conv_cor <- cor(log(working_tra$G), log(last_tra$G))
            }
            else {
                conv_cor <- cor(log(working_tra), log(last_tra), 
                  "complete.obs")
            }
        }
        if (conv_cor >= cor_threshold) {
            converged <- TRUE
        }
    }
    effects <- vector(mode = "list", length = 3)
    all.dim <- c(1, 2, 3)
    miss.dim <- all.dim[-intersect(effect_order, all.dim)]
    miss_effect <- rep.int(1, dim(tra)[[miss.dim]])
    names(miss_effect) <- dimnames(tra)[[miss.dim]]
    prim_effect <- est_effect(working_tra, effect_order[1], mean_type, 
        sparse)
    sec.series <- remove_effect(tra, prim_effect, effect_order[1], 
        form, sparse)
    sec_effect <- est_effect(sec.series, effect_order[2], mean_type, 
        sparse)
    effects <- list(prim_effect, sec_effect, miss_effect)
    effect_order <- order(c(effect_order[1], effect_order[2], 
        miss.dim))
    effects <- effects[effect_order]
    effects <- sort_effects(effects, tra, sparse)
    effects <- rescale_effects(effects, form)
    fit <- model_fit_tra(effects, tra, model, form, error, sparse)
    settings <- list(model = model, form = form, error = error, 
        sparse = sparse, method = "old_sfs")
    out <- list(effects = effects, tra = tra, fit = fit, settings = settings)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
