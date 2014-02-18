\name{pad_effects}
\alias{pad_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pad_effects(effects, tra, form = "multiplicative", sparse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
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
function (effects, tra, form = "multiplicative", sparse = TRUE) 
{
    if (form == "multiplicative") {
        na.value <- 1
    }
    else {
        na.value <- 0
    }
    new_effects <- list(I = NA, T = NA, A = NA)
    if (sparse) {
        for (i in c("I", "T", "A")) {
            j <- c(I = "i", T = "t", A = "a")[i]
            if (length(effects[[i]] > 0)) {
                new_effects[[i]] <- effects[[i]]
                if (length(effects[[i]]) < nlevels(tra[[j]])) {
                  missing_names <- levels(tra[[j]])[!(levels(tra[[j]]) \%in\% 
                    names(effects[[i]]))]
                  missing_effects <- rep(na.value, times = length(missing_names))
                  names(missing_effects) <- missing_names
                  new_effects[[i]] <- c(effects[[i]], missing_effects)
                }
            }
            else {
                tra_dim <- which(c("I", "T", "A") == i)
                new_effects[[i]] <- rep.int(na.value, nlevels(tra[[tra_dim + 
                  1]]))
                effect_names <- levels(tra[[tra_dim + 1]])
                if (i == "T" | i == "A") {
                  effect_names <- as.numeric(as.character(effect_names))
                }
                names(new_effects[[i]]) <- sort(effect_names)
            }
        }
    }
    else {
        for (i in c("I", "T", "A")) {
            j <- c(I = 1, T = 2, A = 3)[i]
            if (length(effects[[i]] > 0)) {
                new_effects[[i]] <- effects[[i]]
                if (length(effects[[i]]) < dim(tra)[j]) {
                  missing_names <- dimnames(tra)[[j]][!(dimnames(tra)[[j]] \%in\% 
                    names(effects[[i]]))]
                  missing_effects <- rep(na.value, times = length(missing_names))
                  names(missing_effects) <- missing_names
                  new_effects[[i]] <- c(effects[[i]], missing_effects)
                }
            }
            else {
                tra_dim <- which(c("I", "T", "A") == i)
                new_effects[[i]] <- rep.int(na.value, dim(tra)[tra_dim])
                names(new_effects[[i]]) <- dimnames(tra)[[tra_dim]]
            }
        }
    }
    return(new_effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
