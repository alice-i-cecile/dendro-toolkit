\name{make_standardization_plots}
\alias{make_standardization_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_standardization_plots(effects, data, link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{
%%     ~~Describe \code{effects} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (effects, data, link = "log") 
{
    plots <- list()
    sample_depth_time_my_plot <- make_sample_depth_plot(data$original, 
        "Time")
    sample_depth_age_my_plot <- make_sample_depth_plot(data$original, 
        "Age")
    plots <- c(plots, list(sample_depth_time_plot = sample_depth_time_plot, 
        sample_depth_age_plot = sample_depth_age_plot))
    if ("Tree" \%in\% names(effects)) {
        tree_effect_my_plot <- make_effect_plot(effects, "Tree")
        tree_effect_density_my_plot <- make_effect_density_plot(effects, 
            "Tree", link)
        tree_effect_age_my_plot <- make_tree_effect_age_plot(effects, 
            data$original)
        tree_effect_year_my_plot <- make_tree_effect_year_plot(effects, 
            data$original)
        plots <- c(plots, list(tree_effect_plot = tree_effect_plot, 
            tree_effect_density_plot = tree_effect_density_plot))
    }
    if ("Time" \%in\% names(effects)) {
        time_effect_my_plot <- make_effect_plot(effects, "Time", 
            temporal = TRUE)
        time_effect_density_my_plot <- make_effect_density_plot(effects, 
            "Time", link)
        plots <- c(plots, list(time_effect_plot = time_effect_plot, 
            time_effect_density_plot = time_effect_density_plot))
    }
    if ("Age" \%in\% names(effects)) {
        age_effect_my_plot <- make_effect_plot(effects, "Age", 
            temporal = TRUE)
        age_effect_density_my_plot <- make_effect_density_plot(effects, 
            "Age", link)
        plots <- c(plots, list(age_effect_plot = age_effect_plot, 
            age_effect_density_plot = age_effect_density_plot))
    }
    residual_density_my_plot <- make_residual_density_plot(data$residuals, 
        link)
    plots <- c(plots, residual_density_plot = residual_density_plot)
    return(plots)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
