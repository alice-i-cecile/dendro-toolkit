\name{make_sample_depth_plot}
\alias{make_sample_depth_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a plot showing the number of samples in the chronology by year or age
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_sample_depth_plot(tra, id = "Time")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
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
function (tra, id = "Time") 
{
    sample_depth <- sample_depth_tra(tra, id)
    df <- data.frame(sample_depth = sample_depth, id = levels(tra[[id]]))
    df$id <- as.numeric(as.character(df$id))
    my_plot <- ggplot(df, aes(x = id, y = sample_depth)) + geom_area() + 
        xlab(id) + ylab("Sample depth") + theme_bw()
    return(my_plot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
