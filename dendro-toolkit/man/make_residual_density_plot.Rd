\name{make_residual_density_plot}
\alias{make_residual_density_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_residual_density_plot(residuals, link)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{residuals}{
%%     ~~Describe \code{residuals} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
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
function (residuals, link) 
{
    df <- data.frame(residuals = residuals$Growth)
    my_plot <- ggplot(df, aes(x = effect)) + theme_bw() + ylab("Density estimate") + 
        xlab("Residuals") + geom_density(colour = "red", fill = "red", 
        alpha = 0.5)
    x_min <- min(min(scaled_effect), 0)
    x_max <- max(scaled_effect)
    x_ticks <- seq(from = x_min, to = x_max, length.out = 100)
    if (link == "identity") {
        pdf_sd <- sd(effect)
        pdf <- dnorm(x_ticks, sd = pdf_sd)
    }
    else if (link == "log") {
        pdf_sd <- sd(log(effect))
        pdf <- dlnorm(x_ticks, sdlog = pdf_sd)
    }
    pdf_data <- data.frame(effect = x_ticks, density = pdf)
    my_plot <- my_plot + geom_area(data = pdf_data, aes(x = effect, 
        y = density), colour = "black", fill = "black", alpha = 0.5)
    return(my_plot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
