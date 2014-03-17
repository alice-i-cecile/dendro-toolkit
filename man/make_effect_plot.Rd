\name{make_effect_plot}
\alias{make_effect_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot coefficients for a tree-ring standardization model
}
\description{
Produce a plot showing the coefficients for an effect in a tree-ring standardization model
}
\usage{
make_effect_plot(effects, effect_name, temporal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{effect_name}{
%%     ~~Describe \code{effect_name} here~~
}
  \item{temporal}{
%%     ~~Describe \code{temporal} here~~
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
function (effects, effect_name, temporal = FALSE) 
{
    df <- data.frame(effect = effects[[effect_name]], id = names(effects[[effect_name]]))
    if (temporal) {
        df$id <- as.numeric(as.character(df$id))
    }
    my_plot <- ggplot(df, aes(x = id, y = effect)) + theme_bw() + 
        ylab(paste(effect_name, "effect")) + xlab(effect_name)
    if (temporal) {
        my_plot <- my_plot + geom_line()
    }
    else {
        my_plot <- my_plot + geom_bar()
    }
    return(my_plot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
