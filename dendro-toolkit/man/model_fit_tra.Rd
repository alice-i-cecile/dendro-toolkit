\name{model_fit_tra}
\alias{model_fit_tra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model_fit_tra(effects, tra, model, form, error, method = "sfs", k = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
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
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (effects, tra, model, form, error, method = "sfs", k = NA) 
{
    fit <- list()
    if (sum(unlist(model)) == 0) {
        fit$predicted <- tra
        if (form == "additive") {
            fit$predicted$G <- 0
        }
        else {
            fit$predicted$G <- 1
        }
        fit$residuals <- tra
    }
    else {
        fit$predicted <- predicted_tra(effects, tra, form)
        fit$residuals <- residuals_tra(tra, fit$predicted, error)
    }
    fit$n <- n_tra(tra)
    if (method == "gam") {
        fit$k <- k
    }
    else {
        fit$k <- k_tra(tra, model)
    }
    fit$sigma <- sigma_tra(fit$residuals, error)
    fit$tss <- tss_tra(tra, error)
    fit$rss <- rss_tra(fit$residuals, error)
    fit$llh <- llh_tra(fit$residuals, error)
    fit$Rsq <- Rsq_tra(fit$rss, fit$tss)
    fit$adj.Rsq <- adj.Rsq_tra(fit$rss, fit$tss, fit$n, fit$k)
    fit$AIC <- AIC_tra(fit$llh, fit$k)
    fit$AICc <- AICc_tra(fit$llh, fit$k, fit$n)
    fit$BIC <- BIC_tra(fit$llh, fit$k, fit$n)
    return(fit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
