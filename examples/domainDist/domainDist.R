#domainDist
library(chroGPS)

# factors domain dist test to assess differential analysis goodness
load('/Volumes/biostats/research/chroGPS/data/Paper/Main/s2_Jan13.RData')
load('/Volumes/biostats/research/chroGPS/data/Paper/Main/bg3_Jan13.RData')

> s2
RangedDataList of length 76
names(76): ASH1-Q4177.S2 BEAF-70.S2 BEAF-HB.S2 Chro(Chriz)BR.S2 Chro(Chriz)WR.S2 CP190-HB.S2 CP190-VC.S2 ... Su(Hw)-VC.S2 Su(var)3-7-Q3448.S2 Su(var)3-9-Q2598.S2 Su(var)3-9.S2 WDS_Q2691.S2 ZW5.S2

> head(s2[[1]][,1:4])
RangedData with 6 rows and 4 value columns across 8 spaces
                    space             ranges |        peak      strand     feature start_position
                 <factor>          <IRanges> | <character> <character> <character>      <numeric>
0004 FBgn0016977    chr2L [ 158853,  159885] |        0004           + FBgn0016977         159034
0012 FBgn0003963    chr2L [ 523753,  526620] |        0012           + FBgn0003963         523467
0016 FBgn0003926    chr2L [1432327, 1433480] |        0016           + FBgn0003926        1432711
0017 FBgn0086758    chr2L [1649515, 1663227] |        0017           + FBgn0086758        1651260
0018 FBgn0086758    chr2L [1663345, 1665545] |        0018           + FBgn0086758        1651260
0019 FBgn0086758    chr2L [1665677, 1672308] |        0019           + FBgn0086758        1651260

# Unify replicates
mnames <- sort(unique(intersect(s2names$Factor,bg3names$Factor)))
s2.repset <- mclapply(mnames,function(i) do.call(rbind,s2[s2names$Factor==i]),mc.cores=12)
bg3.repset <- mclapply(mnames,function(i) do.call(rbind,bg3[bg3names$Factor==i]),mc.cores=12)
names(s2.repset) <- names(bg3.repset) <- mnames
                       
> head(s2names[,-4])
    ExperimentName      Factor             ID      Color Darkcolor
1    ASH1-Q4177.S2        ASH1 modENCODE_2984 lightgreen darkgreen
2       BEAF-70.S2        BEAF  modENCODE_922       grey  darkgrey
3       BEAF-HB.S2        BEAF  modENCODE_274       grey  darkgrey
4 Chro(Chriz)BR.S2 CHRO(CHRIZ)  modENCODE_278 lightgreen darkgreen
5 Chro(Chriz)WR.S2 CHRO(CHRIZ)  modENCODE_279 lightgreen darkgreen
6      CP190-HB.S2       CP190  modENCODE_925       grey  darkgrey
> 

# Generate unified domain names
domains$Domain <- domains$Color
domains <- unique(s2names[s2names$Factor %in% mnames,c('Factor','Color','Domain')])
domains$Domain[domains$Color %in% c('lightgreen','purple')] <- 'Active'
domains$Domain[domains$Color=='lightblue'] <- 'HP1a repression'
domains$Domain[domains$Color=='yellow'] <- 'Polycomb repression'
domains$Domain[domains$Color=='grey'] <- 'Boundaries'
                       
\begin{lstlisting}[backgroundcolor=\color{black!10}]
# Compute distances
d.s2 <- distGPS(RangedDataList(s2.repset),metric='avgdist',mc.cores=12)
d.bg3 <- distGPS(RangedDataList(bg3.repset),metric='avgdist',mc.cores=12)
s2names.repset <- unique(s2names[s2names$Factor %in% mnames,c('Factor','Color','Domain')])
bg3names.repset <- unique(bg3names[bg3names$Factor %in% mnames,c('Factor','Color','Domain')])

# Compute inter-domain distances
dd.s2 <- domainDist(as.matrix(d.s2),gps='factors',domain=domains$Color,type='inter')
dd.bg3 <- domainDist(as.matrix(d.bg3),gps='factors',domain=domains$Color,type='inter')

# Random seed
set.seed(149)

# Alterate s2
s2.alt <- s2.repset
s2.alt[['EZ']] <- RangedData(rbind(as.data.frame(s2.repset[['GAF']])[sample(1:(nrow(s2.repset[['GAF']])/2)),],as.data.frame(s2.repset[['HP1B']])[sample(1:(nrow(s2.repset[['HP1B']])/2)),]))
d.s2.alt <- distGPS(RangedDataList(s2.alt),metric='avgdist',mc.cores=12)

# Plot S2 vs BG3
par(las=1,mar=c(4,8,4,4))
mycors1 <- rev(diag(cor(as.matrix(d.s2),as.matrix(d.bg3))))
barplot(mycors1,horiz=TRUE,xlim=c(0,1),main='S2 / BG3',col=s2names.repset[names(mycors1),'Color'],font=2)
for (i in 1:length(summary(mycors))) abline(v=summary(mycors)[i],col=i,lwd=2,lty=3)

# Plot S2 Altered vs BG3
par(las=1,mar=c(4,8,4,4))
mycors2 <- rev(diag(cor(as.matrix(d.s2.alt),as.matrix(d.bg3))))
barplot(mycors2,horiz=TRUE,xlim=c(0,1),main='S2 Altered / BG3',col=bg3names.repset[names(mycors2),'Color'],font=2)
for (i in 1:length(summary(mycors))) abline(v=summary(mycors)[i],col=i,lwd=2,lty=3)

## Perform Mantel/Permutation test
## Mantel test is based on comparisons of total and per-domain sub-distance matrixes
## Permutation test is exactly what its name indicates, using regioneR to evaluate loss of overlap, it is, however, highly intensive computationally speaking
## ...
