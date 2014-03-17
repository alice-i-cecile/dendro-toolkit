\name{make_tree_effect_age_plot}
\alias{make_tree_effect_age_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a plot showing tree effect vs. age at sampling
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_tree_effect_age_plot(effects, tra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
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
function (effects, tra) 
{
    ages <- sapply(names(effects$Tree), grab_age)
    df <- data.frame(tree = effects$Tree, age = ages)
    my_plot <- ggplot(df, aes(x = age, y = tree)) + geom_point() + 
        geom_smooth() + theme_bw() + ylab("Tree effect") + xlab("Age at sampling")
    return(my_plot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
