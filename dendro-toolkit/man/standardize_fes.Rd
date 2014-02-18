\name{standardize_fes}
\alias{standardize_fes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_fes(tra, incQ = T, incF = T, incA = T, multiplicative = T, corr_m = T, model_type = "lm", span = "auto", degree = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{incQ}{
%%     ~~Describe \code{incQ} here~~
}
  \item{incF}{
%%     ~~Describe \code{incF} here~~
}
  \item{incA}{
%%     ~~Describe \code{incA} here~~
}
  \item{multiplicative}{
%%     ~~Describe \code{multiplicative} here~~
}
  \item{corr_m}{
%%     ~~Describe \code{corr_m} here~~
}
  \item{model_type}{
%%     ~~Describe \code{model_type} here~~
}
  \item{span}{
%%     ~~Describe \code{span} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
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
function (tra, incQ = T, incF = T, incA = T, multiplicative = T, 
    corr_m = T, model_type = "lm", span = "auto", degree = 1, 
    ...) 
{
    tra[which(tra == 0)] <- NA
    if (multiplicative) {
        tra <- log(tra)
    }
    tab.tra <- tra.to.table(tra)
    tab.tra[!is.finite(tab.tra$G), "G"] <- NA
    if (model_type == "lm") {
        growth_formula <- as.formula(make_lm_formula(incQ, incF, 
            incA))
    }
    else if (model_type == "gnm") {
        growth_formula <- as.formula(make_gnm_formula(incQ, incF, 
            incA))
    }
    else if (model_type == "gam") {
        growth_formula <- as.formula(make_gam_formula(incQ, incF, 
            incA))
    }
    print("Model initialized.")
    print(growth_formula)
    if (model_type == "lm") {
        growth_model <- lm(growth_formula, data = tab.tra, ...)
    }
    else if (model_type == "gnm") {
        growth_model <- gnm(growth_formula, data = tab.tra, ...)
    }
    else if (model_type == "gam") {
        if (incA == FALSE) {
            growth_model <- lm(growth_formula, data = tab.tra, 
                ...)
        }
        else {
            tab.tra$A.num <- as.numeric(as.character(tab.tra$A.cv))
            min_span <- (degree + 1)/length(unique(tab.tra$A.num))
            max_span <- 1
            if (span == "auto") {
                print("Automatically selecting span parameter by AIC.")
                optim_span <- function(span) {
                  local_formula <- growth_formula
                  environment(local_formula) <- environment()
                  model <- gam(local_formula, data = tab.tra)
                  model_AIC <- AIC(model)
                  print(c(Span = span, AIC = model_AIC))
                  return(model_AIC)
                }
                span <- optimize(f = optim_span, interval = c(min_span, 
                  max_span))$minimum
                print(paste("Span of", span, "selected."))
            }
            else {
                if (span < min_span) {
                  span <- min_span
                  print(paste("Selected span too low, automatically set to", 
                    min_span))
                }
                if (span > max_span) {
                  span <- max_span
                  print(paste("Selected span too high, automatically set to", 
                    1))
                }
            }
            growth_model <- gam(growth_formula, data = tab.tra, 
                ...)
        }
    }
    model_summ <- summary(growth_model)
    print("Model fit complete.")
    cv <- extract_cv(growth_model, model_type)
    cv <- pad_cv(cv, tra, multiplicative)
    cv <- sort_cv(cv, tra)
    print("Canonical vectors and standard errors of estimate extracted.")
    if (multiplicative) {
        full_cv <- lapply(cv, length) > 0
        cv[full_cv] <- lapply(cv[full_cv], exp)
    }
    cv <- rescaleCV(cv)
    sd.epsilon <- sd(residuals(growth_model), na.rm = T)
    print(paste("Standard deviation of residuals:", signif(sd.epsilon, 
        3)))
    if (model_type == "lm") {
        Rsq <- model_summ$r.squared
        adj.Rsq <- model_summ$adj.r.squared
    }
    loglike <- logLik(growth_model)
    aic <- AIC(growth_model)
    bic <- BIC(growth_model)
    if (model_type == "lm") {
        fit <- c(Rsq, adj.Rsq, loglike, aic, bic)
        names(fit) <- c("R-squared", "Adjusted.R-squared", "Log-Likelihood", 
            "AIC", "BIC")
    }
    else {
        fit <- c(loglike, aic, bic)
        names(fit) <- c("Log-Likelihood", "AIC", "BIC")
    }
    print("Fit computed.")
    print(fit)
    names(cv) <- c("Q", "F", "A")
    output <- list(cv = cv, RSE = sd.epsilon, fit = fit, model = growth_model)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
