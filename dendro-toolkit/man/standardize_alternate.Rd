\name{standardize_alternate}
\alias{standardize_alternate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
standardize_alternate(tra, model = c("Time", "Age"), link = "log", cor_threshold = 0.999999)
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
function (tra, model = c("Time", "Age"), link = "log", cor_threshold = 0.999999) 
{
    if (link == "log") {
        error <- "lnorm"
        form <- "multiplicative"
    }
    else if (link == "identity") {
        error <- "norm"
        form <- "additive"
    }
    if (error == "lnorm") {
        mean_type <- "geometric"
    }
    else {
        mean_type <- "arithmetic"
    }
    effects <- vector(mode = "list", length = length(model))
    names(effects) <- model
    for (i in model) {
        dim_i <- nlevels(tra[[i]])
        if (form == "additive") {
            effects[[i]] <- rep.int(0, dim_i)
        }
        else {
            effects[[i]] <- rep.int(1, dim_i)
        }
        names(effects[[i]]) <- levels(tra[[i]])
    }
    if (length(model) == 0) {
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
        for (id in model) {
            est_j <- est_effect(working_tra, id, mean_type)
            working_tra <- remove_effect(working_tra, est_j, 
                id, link)
            if (form == "additive") {
                effects[[id]] <- effects[[id]] + est_j
            }
            else {
                effects[[id]] <- effects[[id]] * est_j
            }
        }
        if (error == "norm") {
            conv_cor <- cor(working_tra$Growth, last_tra$Growth)
        }
        else {
            conv_cor <- cor(log(working_tra$Growth), log(last_tra$Growth))
        }
        if (conv_cor >= cor_threshold) {
            converged <- TRUE
        }
    }
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
