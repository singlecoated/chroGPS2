\name{diffGenes}
\alias{diffGenes}
\alias{find.fdr}
\alias{find.threshold}
\alias{diffGenes-methods}
\title{
Performs differential analysis of chroGPS-genes maps based on Bayesian
Correct Classification Rates.
}
\description{
The function uses Bayesian posterior probability of correct
classification (CCR) of chroGPS-genes maps as computed by the
\code{clusGPS} function to identify statistically significant epigenetic
changes between two different backgrounds (technical and/or biological).
The functions \code{find.fdr} and \code{find.threshold} are used
internally to compute approximated FDRs from a vector of input posterior
probabilities of correct classification (CCRs).
}
\usage{
diffGenes(xy,m,clus,label.x,label.y,clusName=NULL,fdr=TRUE,mc.cores=1)
}
\arguments{
  \item{xy}{
    A \code{matrix} or \code{data.frame} object, containing a combined
    matrix with whole-genome epigenetic information for backgrounds X
    and Y, as obtained by the \code{combineGenesMatrix} function.
  }
  \item{m}{
    A \code{mds} object, containing a chroGPS-genes MDS map, obtained
    for the combined matrix xy.
  }
  \item{clus}{
    A \code{clusGPS} object, containing a chroGPS-genes cluster solution, obtained
    for the provided background via the \code{clusGPS} and
    \code{mergeClusters} functions.
  }
  \item{label.x}{
    Name for background X, by default 'mds1'
  }
  \item{label.y}{
    Name for background Y, by default 'mds2'
  }
  \item{clusName}{
    If \code{clus} contains more than one clustering solution
    (i.e. \code{cutree} dendrogram height cut), the name of the selected one. If no name is provided, the first element is used.
  }
  \item{fdr}{
    Set to TRUE (default) to compute False Discovery Rate of cluster
    correct classification from posterior probabilities.
  }
  \item{mc.cores}{
    Number of cores to use for parallel computations.
  }
}
\value{
The function returns a data frame with posterior probabilities and FDRs
of cluster classifications for all genes in backgrounds X and Y, so that
significant changes can be traced and reported easily.
}
\author{
Oscar Reina.
}
\seealso{
See functions \code{distGPS, mds, clusGPS, profileClusters} for
generating chroGPS-genes maps. See function \code{plotTransitions} for
examples on how to visualize differential gene maps results.
}
\examples{
}
\keyword{ ~mds }
\keyword{ ~clustering }
