\name{standardize_tra}
\alias{standardize_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_tra(tra, model = list(I = FALSE, T = TRUE, A = TRUE), form = "multiplicative", error = "lnorm", method = "sfs", post_hoc = TRUE, ...)
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
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{post_hoc}{
%%     ~~Describe \code{post_hoc} here~~
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
    error = "lnorm", method = "sfs", post_hoc = TRUE, ...) 
{
    if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra <= 
        0], na.rm = TRUE) > 0)) {
        if (form == "multiplicative") {
            warning("Zero or negative values found in data when using a multiplicative model. Estimated effects will not be stable.")
        }
        if (error == "lnorm") {
            warning("Zero or negative values found in data when using a lognormal error form. Model fitting will fail, as non-positive values can not be generated using a log-normal distribution.")
        }
    }
    if ((form == "additive" & error == "lnorm") | (form == "multiplicative" & 
        error == "norm")) {
        warning("Model form and error distribution are an unrealistic match. Be sure to check the model residuals. Model fitting problems may arise.")
    }
    if (method == "likelihood") {
        out <- standardize_likelihood(tra, model, form, error, 
            ...)
    }
    else if (method == "least_squares") {
        out <- standardize_least_squares(tra, model, form, error, 
            ...)
    }
    else if (method == "sfs") {
        out <- standardize_tsfs(tra, model, form, error, ...)
    }
    else if (method == "rcs") {
        out <- standardize_rcs(tra, model, form, error, ...)
    }
    else if (method == "gam") {
        out <- standardize_gam(tra, model, form, error, ...)
    }
    if (sum(unlist(model)) == 3) {
        if (post_hoc) {
            out$effects <- post_hoc_intercession(out$effects, 
                out$tra, form)
            warning("Three effect model selected. Post-hoc selection was used to stabilize parameter estimates.")
        }
        else {
            warning("Three effect model selected. Parameter estimates are wildly unreliable. Consider using post-hoc selection.")
        }
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
