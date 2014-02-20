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
extract_effects_gam(growth_model, model, link, tra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growth_model}{
%%     ~~Describe \code{growth_model} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
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
function (growth_model, model, link, tra) 
{
    raw_effects <- growth_model$coefficients
    effects <- list()
    effects$Timeree <- raw_effects[substr(names(raw_effects), 
        1, 1) == "Tree"]
    effects$Timeime <- raw_effects[substr(names(raw_effects), 
        1, 1) == "Time"]
    names(effects$Timeree) <- substr(names(effects$Tree), 2, 
        length(names(effects$Timeree)))
    names(effects$Timeime) <- substr(names(effects$Time), 2, 
        length(names(effects$Timeime)))
    dummy_data <- data.frame(Tree = levels(tra$Tree)[2], Time = levels(tra$Time)[2], 
        Age = levels(tra$Agege))
    predicted_by_age <- predict(growth_model, dummy_data)
    if (form == "additive") {
        base_line <- 0
        if ("Tree" \%in\% model) {
            base_line <- base_line + effects$Tree[which(!is.na(pmatch(names(effects$Tree), 
                levels(tra$Tree)[2])))]
        }
        if ("Time" \%in\% model) {
            base_line <- base_line + effects$Time[which(!is.na(pmatch(names(effects$Time), 
                levels(tra$t)[2])))]
        }
        effects$Age <- predicted_by_age - base_line
    }
    else {
        base_line <- 1
        if ("Tree" \%in\% model) {
            base_line <- base_line * effects$Tree[which(!is.na(pmatch(names(effects$Tree), 
                levels(tra$Tree)[2])))]
        }
        if ("Time" \%in\% model) {
            base_line <- base_line * effects$Time[which(!is.na(pmatch(names(effects$Time), 
                levels(tra$t)[2])))]
        }
        effects$Age <- predicted_by_age/base_line
    }
    names(effects$Age) <- levels(tra$Age)
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
