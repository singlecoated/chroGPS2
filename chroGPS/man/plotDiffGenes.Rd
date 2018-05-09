\name{plotDiffGenes}
\alias{plotDiffGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the result(s) from the differential analysis of chroGPS-genes maps.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDiffGenes(x, m, clus, res, label.x = "1", label.y = "2", clusName = NULL, transitions = NULL, plotIDs = TRUE, fdr1 = 0.05, fdr2 = 0.05, minGenes = 1, plotCentroids = TRUE, centroid.lwd = 2, centroid.lty = 1, point.col = "lightgrey", probContour = 0.75, contour.lwd = 2, contour.lty = 2, dolegend = TRUE, poslegend = "bottomright", xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{clus}{
%%     ~~Describe \code{clus} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{label.x}{
%%     ~~Describe \code{label.x} here~~
}
  \item{label.y}{
%%     ~~Describe \code{label.y} here~~
}
  \item{clusName}{
%%     ~~Describe \code{clusName} here~~
}
  \item{transitions}{
%%     ~~Describe \code{transitions} here~~
}
  \item{plotIDs}{
%%     ~~Describe \code{plotIDs} here~~
}
  \item{fdr1}{
%%     ~~Describe \code{fdr1} here~~
}
  \item{fdr2}{
%%     ~~Describe \code{fdr2} here~~
}
  \item{minGenes}{
%%     ~~Describe \code{minGenes} here~~
}
  \item{plotCentroids}{
%%     ~~Describe \code{plotCentroids} here~~
}
  \item{centroid.lwd}{
%%     ~~Describe \code{centroid.lwd} here~~
}
  \item{centroid.lty}{
%%     ~~Describe \code{centroid.lty} here~~
}
  \item{point.col}{
%%     ~~Describe \code{point.col} here~~
}
  \item{probContour}{
%%     ~~Describe \code{probContour} here~~
}
  \item{contour.lwd}{
%%     ~~Describe \code{contour.lwd} here~~
}
  \item{contour.lty}{
%%     ~~Describe \code{contour.lty} here~~
}
  \item{dolegend}{
%%     ~~Describe \code{dolegend} here~~
}
  \item{poslegend}{
%%     ~~Describe \code{poslegend} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
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
function (x, m, clus, res, label.x = "1", label.y = "2", clusName = NULL, 
    transitions = NULL, plotIDs = TRUE, fdr1 = 0.05, fdr2 = 0.05, 
    minGenes = 1, plotCentroids = TRUE, centroid.lwd = 2, centroid.lty = 1, 
    point.col = "lightgrey", probContour = 0.75, contour.lwd = 2, 
    contour.lty = 2, dolegend = TRUE, poslegend = "bottomright", 
    xlim, ylim, ...) 
{
    if (is.null(clusName)) 
        clusName <- clusNames(clus)[1]
    if (missing(xlim)) 
        xlim <- range(m@points)
    if (missing(ylim)) 
        ylim <- range(m@points)
    fdr <- res[, grep("FDR", colnames(res))]
    if (!ncol(fdr) \%in\% c(0, 2)) 
        stop("Invalid Differential results table")
    if (ncol(fdr) == 2) 
        res <- res[fdr[, 1] < fdr1 & fdr[, 2] < fdr2, ]
    res$cc <- apply(res[, grep("ClusID", colnames(res))], 1, 
        paste, collapse = ".")
    if (is.null(transitions)) 
        transitions <- sort(unique(res$cc))
    for (mycc in transitions) {
        plot(m, point.col = point.col, xlim = xlim, ylim = ylim, 
            ...)
        plot(clus, type = "contours", k = clusName, lwd = contour.lwd, 
            probContour = probContour, lty = contour.lty, xlim = xlim, 
            ylim = ylim, ...)
        if (sum(res$cc == mycc) >= minGenes) 
            plotTransitions(x, m, clus, res, plotCentroids = TRUE, 
                transID = mycc, xlim = xlim, ylim = ylim)
        if (dolegend) {
            tclus <- tabClusters(clus, clusName)
            legend(poslegend, legend = paste(names(tclus), as.numeric(tclus), 
                sep = ":"), col = rainbow(length(tclus)), pch = 19, 
                bty = "n", cex = 1, title = "Cluster: # of elements")
        }
        plotCenters(m, clus, clusName = clusName)
        legend("topleft", legend = sprintf("r2=\%.3f / stress=\%.3f", 
            m@R.square, m@stress))
        legend("topright", legend = sprintf("Cluster Transition: \%s\nFDR=\%.2f/\%.2f (n=\%d)", 
            mycc, fdr1, fdr2, sum(res$cc == mycc)))
        legend("bottomleft", col = c("blue", "red"), lwd = 2, 
            legend = c(label.x, label.y))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
