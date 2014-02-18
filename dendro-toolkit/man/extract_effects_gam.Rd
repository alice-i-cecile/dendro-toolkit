\name{extract_effects_gam}
\alias{extract_effects_gam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_effects_gam(growth_model, model, form, tra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growth_model}{
%%     ~~Describe \code{growth_model} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
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
function (growth_model, model, form, tra) 
{
    raw_effects <- growth_model$coefficients
    effects <- list()
    effects$I <- raw_effects[substr(names(raw_effects), 1, 1) == 
        "i"]
    effects$T <- raw_effects[substr(names(raw_effects), 1, 1) == 
        "t"]
    names(effects$I) <- substr(names(effects$I), 2, length(names(effects$I)))
    names(effects$T) <- substr(names(effects$T), 2, length(names(effects$T)))
    dummy_data <- data.frame(i = levels(tra$i)[2], t = levels(tra$t)[2], 
        a = levels(tra$a))
    predicted_by_age <- predict(growth_model, dummy_data)
    if (form == "additive") {
        base_line <- 0
        if (model$I) {
            base_line <- base_line + effects$I[which(!is.na(pmatch(names(effects$I), 
                levels(tra$i)[2])))]
        }
        if (model$T) {
            base_line <- base_line + effects$T[which(!is.na(pmatch(names(effects$T), 
                levels(tra$t)[2])))]
        }
        effects$A <- predicted_by_age - base_line
    }
    else {
        base_line <- 1
        if (model$I) {
            base_line <- base_line * effects$I[which(!is.na(pmatch(names(effects$I), 
                levels(tra$i)[2])))]
        }
        if (model$T) {
            base_line <- base_line * effects$T[which(!is.na(pmatch(names(effects$T), 
                levels(tra$t)[2])))]
        }
        effects$A <- predicted_by_age/base_line
    }
    names(effects$A) <- levels(tra$a)
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
