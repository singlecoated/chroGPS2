\name{rankFactors}
\alias{rankFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to help selecting candidate epigenetic factors based on
Epigenetic Domain cohesion/separation or global predictive properties.  
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankFactors(d, sampleinfo, ranktype = "domainDist", selName = "Color", selValue, k = NULL, mc.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{sampleinfo}{
%%     ~~Describe \code{sampleinfo} here~~
}
  \item{ranktype}{
%%     ~~Describe \code{ranktype} here~~
}
  \item{selName}{
%%     ~~Describe \code{selName} here~~
}
  \item{selValue}{
%%     ~~Describe \code{selValue} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
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
function (d, sampleinfo, ranktype = "domainDist", selName = "Color", 
    selValue, k = NULL, mc.cores = 1) 
{
    domains <- sampleinfo
    if (!selName \%in\% colnames(domains)) 
        stop("Invalid sel name in sampleinfo")
    ff <- rownames(domains)[domains$Domain == selValue]
    ff <- combs(ff, k = k)
    cat(sprintf("\n### Evaluating \%s domain ###\n", selValue))
    dd <- parallel::mclapply(1:nrow(ff), function(i) {
        sel <- rownames(domains) \%in\% ff[i, ] | domains$Domain != 
            selValue
        inter <- domainDist(as.matrix(d)[sel, sel], gps = "factors", 
            domain = domains$Domain[sel], type = "inter", plot = FALSE, 
            avg = FALSE)
        inter <- inter[grep(selValue, names(inter))]
        inter <- mean(unlist(inter))
        intra <- mean(domainDist(as.matrix(d)[sel, sel], gps = "factors", 
            domain = domains$Domain[sel], type = "intra", plot = FALSE, 
            avg = TRUE)$Avg)
        return(list(inter = inter, intra = intra))
    }, mc.cores = mc.cores)
    names(dd) <- apply(ff, 1, paste, collapse = " + ")
    dd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
