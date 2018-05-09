\name{profileClusters}
\alias{profileClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assess epigenetic profiles for genes present in each
cluster as obtained by the \code{clusGPS} function.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
profileClusters(x, clus, clusName = NULL, normalize = FALSE, mc.cores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{clus}{
%%     ~~Describe \code{clus} here~~
}
  \item{clusName}{
%%     ~~Describe \code{clusName} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, clus, clusName = NULL, normalize = FALSE, mc.cores = 1, 
    ...) 
{
    x <- as.data.frame(x)
    uid <- NULL
    txt <- paste("uid <- paste(", paste("x[,", 1:ncol(x), "]", 
        collapse = ","), ", sep=',')", sep = "")
    eval(parse(text = txt))
    if (is.null(clusName)) 
        clusName <- clusNames(clus)[1]
    clusID <- clus@clus[[clusName]]$id[uid]
    if (normalize) 
        ans <- do.call(rbind, parallel::mclapply(sort(unique(clusID)), 
            function(i) {
                x.sel <- x[clusID == i, ]
                log2(1 + colMeans(x.sel)) - log2(1 + colMeans(x))
            }, mc.cores = mc.cores))
    else ans <- do.call(rbind, parallel::mclapply(sort(unique(clusID)), 
        function(i) {
            x.sel <- x[clusID == i, ]
            log2(1 + (colMeans(x.sel)))
        }, mc.cores = mc.cores))
    rownames(ans) <- sort(unique(clusID))
    return(ans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
