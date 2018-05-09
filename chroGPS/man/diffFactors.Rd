\name{diffFactors}
\alias{diffFactors}
\alias{diffFactors-methods}
\title{
Performs differential analysis of chroGPS-factors maps based on
Procrustes analysis.
}
\description{
The function uses Procrustes analysis to compare two chroGPS-factors MDS
solutions, providing a visual map highlighting differences and a ranked
list of Procrustes squared errors between analog replicated factors in
both maps.
}
\usage{
diffFactors(m1,m2,name1='mds1',name2='mds2',minPoints=10,poslegend='topleft',plot=TRUE,pointcol.m1='red',pointcol.m2='green',textcol.m1='black',textcol.m2='grey',pch.m1=19,pch.m2=19,cex.m1=1,cex.m2=1,segcol='blue',...)
}
\arguments{
  \item{m1}{
    A \code{mds} object, containing a chroGPS-factors MDS map.
  }
  \item{m2}{
    A \code{mds} object, containing a chroGPS-factors MDS map, must be
    conformable with m1 (same number of elements and in the same order in
    the MDS points).
  }
  \item{name1}{
    Name for m1, by default 'mds1'
  }
  \item{name2}{
    Name for m2, by default 'mds2'
  }
  \item{minPoints}{
    Minimum number of points in each of the MDS objects.
  }
  \item{poslegend}{
    Position for the map legend, passed on to \code{legend} function.
  }
  \item{plot}{
    Wether to plot the differential map or not (in case we only want to
    retrieve Procrustes errors).
  }
  \item{pointcol.m1}{
    Color vector for points in MDS1.
  }
  \item{pointcol.m2}{
    Color vector for points in MDS2.
  }
  \item{textcol.m1}{
    Color vector for labels in MDS1. Labels are taken from MDS points rownames.
  }
  \item{textcol.m2}{
    Color vector for labels in MDS2. Labels are taken from MDS points rownames.
  }
  \item{pch.m1}{
    Pch for MDS1 points.
  }
  \item{pch.m2}{
    Pch for MDS2 points.
  }
  \item{cex.m1}{
    Cex for MDS1 points.
  }
  \item{cex.m2}{
    Cex for MDS2 points.
  }
  \item{segcol}{
    Color for segments joining replicates from both MDS.
  }
  \item{...}{
    Further arguments passed on to the \code{plot} function.
  }
}
\value{
The function returns a data frame with Procrustes errors between paired
replicates in both chroGPS-factors MDS maps.
}
\author{
Oscar Reina.
}
\seealso{
See functions \code{distGPS, mds} for generating chroGPS-factors maps.
}
\examples{
}
\keyword{ ~mds }
