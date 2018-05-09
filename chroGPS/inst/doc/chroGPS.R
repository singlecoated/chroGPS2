### R code from vignette source 'chroGPS.Rnw'

# Code to reprocess data as needed
library(chroGPS)
library(BSgenome.Dmelanogaster.UCSC.dm3)

# s2
load('/Volumes/biostats/routines/R/chroGPS_library/chroGPS/data/s2.RData')
s2 <- lapply(s2@listData,as,'GRanges')
s2 <- lapply(s2,function(x) { seqlengths(x) <- seqlengths(Dmelanogaster)[names(seqlengths(x))]; x })
s2 <- GRangesList(s2)
save(s2,s2names,s2.tab,file='../data/s2.RData')

# bg3
load('/Volumes/biostats/research/chroGPS/data/Paper/Main/bg3_Jan13.RData')
bg3.tab <- bg3.tab[rownames(bg3.tab) %in% rownames(s2.tab),]
sel <- bg3names$Factor %in% s2names$Factor
bg3names <- bg3names[sel,]
bg3 <- lapply(bg3@listData,as,'GRanges')
bg3 <- bg3[sel]
bg3 <- lapply(bg3,function(x) { seqlengths(x) <- seqlengths(Dmelanogaster)[names(seqlengths(x))]; x })
bg3 <- GRangesList(bg3)
save(bg3,bg3names,file='../data/bg3.RData')

# s2Seq
load('/Volumes/biostats/routines/R/chroGPS_library/chroGPS/data/s2Seq.RData')
s2Seq <- lapply(s2Seq@listData,as,'GRanges')
s2Seq <- lapply(s2Seq,function(x) { seqlengths(x) <- seqlengths(Dmelanogaster)[names(seqlengths(x))]; x })
s2Seq <- GRangesList(s2Seq)
save(s2Seq,s2SeqNames,file='../data/s2Seq.RData')

###################################################
### code chunk number 1: import1
###################################################
options(width=70)
par(mar=c(2,2,2,2))
library(chroGPS)
data(s2) # Loading Dmelanogaster S2 modEncode toy example
data(s2Seq)
data(bg3)
data(toydists) # Loading precomputed distGPS objects

###################################################
### code chunk number 2: mds1
###################################################
source('../R/distGPS.R')
d <- distGPS(s2, metric='avgdist')
d
mds1 <- mds(d,k=2,type='isoMDS')
mds1
mds1.3d <- mds(d,k=3,type='isoMDS')
mds1.3d

###################################################
### code chunk number 3: figmds1
###################################################
cols <- as.character(s2names$Color)
plot(mds1,drawlabels=TRUE,point.pch=20,point.cex=8,text.cex=.7,
point.col=cols,text.col='black',labels=s2names$Factor,font=2)
legend('topleft',legend=sprintf('R2=%.3f / stress=%.3f',getR2(mds1),getStress(mds1)),
bty='n',cex=1)
#plot(mds1.3d,drawlabels=TRUE,type.3d='s',point.pch=20,point.cex=.1,text.cex=.7,
#point.col=cols,text.col='black',labels=s2names$Factor)

###################################################
### code chunk number 4: figmds1
###################################################
cols <- as.character(s2names$Color)
plot(mds1,drawlabels=TRUE,point.pch=20,point.cex=8,text.cex=.7,
point.col=cols,text.col='black',labels=s2names$Factor,font=2)
legend('topleft',legend=sprintf('R2=%.3f / stress=%.3f',getR2(mds1),getStress(mds1)),
bty='n',cex=1)
#plot(mds1.3d,drawlabels=TRUE,type.3d='s',point.pch=20,point.cex=.1,text.cex=.7,
#point.col=cols,text.col='black',labels=s2names$Factor)

###################################################
### code chunk number 5: figprocrustes1
###################################################
s2Seq
s2.all <- suppressWarnings(GRangesList(c(lapply(s2,reduce),lapply(s2Seq,reduce)))) # do not complain about chr4
names(s2.all)
d2 <- distGPS(s2.all,metric='avgdist')
mds2 <- mds(d2,k=2,type='isoMDS')
cols <- c(as.character(s2names$Color),as.character(s2SeqNames$Color))
sampleid <- c(as.character(s2names$Factor),as.character(s2SeqNames$Factor))
pchs <- rep(c(20,17),c(length(s2),length(s2Seq)))
point.cex <- rep(c(8,5),c(length(s2),length(s2Seq)))
par(mar=c(2,2,2,2))
plot(mds2,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
legend('topleft',legend=sprintf('R2=%.3f / stress=%.3f',getR2(mds2),getStress(mds2)),
bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

###################################################
### code chunk number 7: figprocrustes2
###################################################
adjust <- rep(c('chip','seq'),c(length(s2),length(s2Seq)))
sampleid <- c(as.character(s2names$Factor),as.character(s2SeqNames$Factor))
mds3 <- procrustesAdj(mds2,d2,adjust=adjust,sampleid=sampleid)
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds3,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
legend('topleft',legend=sprintf('R2=%.3f / stress=%.3f',getR2(mds3),getStress(mds3)),
bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

###################################################
### code chunk number 8: figpeakwidth1
###################################################
s2.pAdj <- suppressWarnings(adjustPeaks(s2.all,adjust=adjust,sampleid=sampleid,logscale=TRUE)) # chr4
d3 <- distGPS(s2.pAdj,metric='avgdist')
mds4 <- mds(d3,k=2,type='isoMDS')
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds4,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
legend('topleft',legend=sprintf('R2=%.3f / s=%.3f',getR2(mds4),getStress(mds4)),
bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

###################################################
### code chunk number 11: import2
###################################################

s2.tab[1:10,1:4]

###################################################
### code chunk number 12: mds6
###################################################
d <- distGPS(s2.tab, metric='tanimoto', uniqueRows=TRUE)
d
mds1 <- mds(d,k=2,type='isoMDS')
mds1
mds2 <- mds(d,k=3,type='isoMDS')
mds2

###################################################
### code chunk number 13: figmds6
###################################################
par(mar=c(2,2,2,2))
plot(mds1,point.cex=1.5,point.col=densCols(getPoints(mds1)))
#plot(mds2,point.cex=1.5,type.3d='s',point.col=densCols(getPoints(mds2)))

###################################################
### code chunk number 14: uniqueCount
###################################################
dim(s2.tab)
dim(uniqueCount(s2.tab))


###################################################
### code chunk number 15: genomeGPS1
###################################################
system.time(mds3 <- mds(d,k=2,type='isoMDS'))
mds3


###################################################
### code chunk number 16: genomeGPS2
###################################################
system.time(mds3 <- mds(d,type='isoMDS',splitMDS=TRUE,split=.5,overlap=.05,mc.cores=4))
mds3
system.time(mds4 <- mds(d,mds3,type='boostMDS',scale=TRUE))
mds4


###################################################
### code chunk number 17: getxpr
###################################################
summary(s2.wt$epigene)
summary(s2.wt$gene)


###################################################
### code chunk number 18: figmds7
###################################################
plot(mds1,point.cex=1.5,scalecol=TRUE,scale=s2.wt$epigene,
     palette=rev(heat.colors(100)))


###################################################
### code chunk number 19: figmds7
###################################################
plot(mds1,point.cex=1.5,scalecol=TRUE,scale=s2.wt$epigene,
     palette=rev(heat.colors(100)))


###################################################
### code chunk number 20: figclus00
###################################################
h <- hclust(as.dist(as.matrix(d)),method='average')
set.seed(149) # Random seed for the MCMC process within density estimation
clus <- clusGPS(d,mds1,h,ngrid=1000,densgrid=FALSE,verbose=TRUE,
preMerge=TRUE,k=max(cutree(h,h=0.5)),minpoints=20,mc.cores=1)
clus


###################################################
### code chunk number 21: figclus0
###################################################
clus
clusNames(clus)
tabClusters(clus,125)
point.col <- rainbow(length(tabClusters(clus,125)))
names(point.col) <- names(tabClusters(clus,125))
point.col
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.col=point.col[as.character(clusterID(clus,125))],
point.pch=19)


###################################################
### code chunk number 22: figclus1
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50))
plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)


###################################################
### code chunk number 23: figclus0
###################################################
clus
clusNames(clus)
tabClusters(clus,125)
point.col <- rainbow(length(tabClusters(clus,125)))
names(point.col) <- names(tabClusters(clus,125))
point.col
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.col=point.col[as.character(clusterID(clus,125))],
point.pch=19)


###################################################
### code chunk number 24: figclus1
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50))
plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)


###################################################
### code chunk number 25: figclus2
###################################################
plot(clus,type='stats',k=max(cutree(h,h=0.5)),ylim=c(0,1),col=point.col,cex=2,pch=19,
lwd=2,ylab='CCR',xlab='Cluster ID',cut=0.75,cut.lty=3,axes=FALSE)
axis(1,at=1:length(tabClusters(clus,125)),labels=names(tabClusters(clus,125))); axis(2)
box()


###################################################
### code chunk number 26: figclus2
###################################################
plot(clus,type='stats',k=max(cutree(h,h=0.5)),ylim=c(0,1),col=point.col,cex=2,pch=19,
lwd=2,ylab='CCR',xlab='Cluster ID',cut=0.75,cut.lty=3,axes=FALSE)
axis(1,at=1:length(tabClusters(clus,125)),labels=names(tabClusters(clus,125))); axis(2)
box()


###################################################
### code chunk number 27: figloc1
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.5,0.95)) plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)
fgenes <- uniqueCount(s2.tab)[,'HP1a_wa184.S2']==1
set.seed(149)
c1 <- contour2dDP(getPoints(mds1)[fgenes,],ngrid=1000,contour.type='none')
for (p in seq(0.1,0.9,0.1)) plotContour(c1,probContour=p,col='black')
legend('topleft',lwd=1,lty=1,col='black',legend='HP1a contours (10 to 90 percent)',bty='n')


###################################################
### code chunk number 28: figloc1
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.5,0.95)) plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)
fgenes <- uniqueCount(s2.tab)[,'HP1a_wa184.S2']==1
set.seed(149)
c1 <- contour2dDP(getPoints(mds1)[fgenes,],ngrid=1000,contour.type='none')
for (p in seq(0.1,0.9,0.1)) plotContour(c1,probContour=p,col='black')
legend('topleft',lwd=1,lty=1,col='black',legend='HP1a contours (10 to 90 percent)',bty='n')


###################################################
### code chunk number 29: figloc2
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.5,0.95)) plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)
set.seed(149) # Random seed for random gene sampling
geneset <- sample(rownames(s2.tab),10,rep=FALSE)
mds2 <- geneSetGPS(s2.tab,mds1,geneset,uniqueCount=TRUE)
points(getPoints(mds2),col='black',cex=5,lwd=4,pch=20)
points(getPoints(mds2),col='white',cex=4,lwd=4,pch=20)
text(getPoints(mds2)[,1],getPoints(mds2)[,2],1:nrow(getPoints(mds2)),cex=1.5)
legend('bottomright',col='black',legend=paste(1:nrow(getPoints(mds2)),
geneset,sep=': '),cex=1,bty='n')


###################################################
### code chunk number 30: figloc2
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.5,0.95)) plot(clus,type='contours',k=max(cutree(h,h=0.5)),lwd=5,probContour=p,
drawlabels=TRUE,labcex=2,font=2)
set.seed(149) # Random seed for random gene sampling
geneset <- sample(rownames(s2.tab),10,rep=FALSE)
mds2 <- geneSetGPS(s2.tab,mds1,geneset,uniqueCount=TRUE)
points(getPoints(mds2),col='black',cex=5,lwd=4,pch=20)
points(getPoints(mds2),col='white',cex=4,lwd=4,pch=20)
text(getPoints(mds2)[,1],getPoints(mds2)[,2],1:nrow(getPoints(mds2)),cex=1.5)
legend('bottomright',col='black',legend=paste(1:nrow(getPoints(mds2)),
geneset,sep=': '),cex=1,bty='n')


###################################################
### code chunk number 31: figmerge0
###################################################
set.seed(149) # Random seed for MCMC within the density estimate process
clus2 <- clusGPS(d,mds1,h,ngrid=1000,densgrid=FALSE,verbose=TRUE,
preMerge=TRUE,k=max(cutree(h,h=0.2)),minpoints=20,mc.cores=1)


###################################################
### code chunk number 32: figmerge1
###################################################
par(mar=c(2,2,2,2))
clus3 <- mergeClusters(clus2,brake=0,mc.cores=1)
clus3
tabClusters(clus3,330)


###################################################
### code chunk number 33: figmerge1
###################################################
par(mar=c(2,2,2,2))
clus3 <- mergeClusters(clus2,brake=0,mc.cores=1)
clus3
tabClusters(clus3,330)

###################################################
### code chunk number 34: figmerge21
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50)) plot(clus2,type='contours',k=max(cutree(h,h=0.2)),
lwd=5,probContour=p,drawlabels=TRUE,labcex=2,font=2)


###################################################
### code chunk number 35: figmerge22
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50)) plot(clus3,type='contours',k=max(cutree(h,h=0.2)),
lwd=5,probContour=p)


###################################################
### code chunk number 36: figmerge21
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50)) plot(clus2,type='contours',k=max(cutree(h,h=0.2)),
lwd=5,probContour=p,drawlabels=TRUE,labcex=2,font=2)


###################################################
### code chunk number 37: figmerge22
###################################################
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds1,point.cex=1.5,point.col='grey')
for (p in c(0.95, 0.50)) plot(clus3,type='contours',k=max(cutree(h,h=0.2)),
lwd=5,probContour=p)


###################################################
### code chunk number 38: figmerge31
###################################################
plot(clus2,type='stats',k=max(cutree(h,h=0.2)),ylim=c(0,1),lwd=2,
ylab='CCR',xlab='Cluster ID')


###################################################
### code chunk number 39: figmerge32
###################################################
plot(clus3,type='stats',k=max(cutree(h,h=0.2)),ylim=c(0,1),lwd=2,
ylab='CCR',xlab='Cluster ID')


###################################################
### code chunk number 40: figmerge31
###################################################
plot(clus2,type='stats',k=max(cutree(h,h=0.2)),ylim=c(0,1),lwd=2,
ylab='CCR',xlab='Cluster ID')


###################################################
### code chunk number 41: figmerge32
###################################################
plot(clus3,type='stats',k=max(cutree(h,h=0.2)),ylim=c(0,1),lwd=2,
ylab='CCR',xlab='Cluster ID')


###################################################
### code chunk number 42: figprofile1
###################################################
source('../R/profileClusters.R')
#p1.1 <- profileClusters(s2.tab, uniqueCount = TRUE, clus=clus3, i=max(cutree(h,h=0.2)), log2 = TRUE, plt = FALSE, minpoints=0)
p1 <- profileClusters2(s2.tab, clus=clus3, clusName=clusNames(clus3),normalize=TRUE)
# Requires gplots library
library(pheatmap)
#pheatmap(log2(p1.1[,1:20]+1),trace='none',col=bluered(100),margins=c(10,12),symbreaks=TRUE,Rowv=FALSE,Colv=FALSE,dendrogram='none')
pheatmap(p1[,1:20],trace='none',col=bluered(100),margins=c(10,12))

###################################################
### code chunk number 44: cyto1
###################################################
# For instance if mds1 contains a valid chroGPS-factors map.
# gps2xgmml(mds1, fname='chroGPS_factors.xgmml', fontSize=4,
# col=s2names$Color, cex=8)
# And use Cytoscape -> File -> Import -> Network (Multiple File Types) 
# to load the generated .xgmml file

###################################################
### code chunk number 45: getData
### https://github.com/singlecoated/chroGPS2/tree/master/examples/getData
###################################################

###################################################
### code chunk number 46: domainDist
### https://github.com/singlecoated/chroGPS2/tree/master/examples/domainDist
###################################################

# Unify replicates
source('../R/mergeReplicates.R')
mnames <- sort(unique(intersect(s2names$Factor,bg3names$Factor)))
sel <- s2names$Factor %in% mnames
s2.repset <- mergeReplicates(s2[sel],id=s2names$Factor[sel],mergeBy='any')
sel <- bg3names$Factor %in% mnames
bg3.repset <- mergeReplicates(bg3[sel],id=bg3names$Factor[sel],mergeBy='any')
# Generate unified domain names
color2domain <- c('Active','Active','HP1a','Polycomb','Boundaries')
names(color2domain) <- c('lightgreen','purple','lightblue','yellow','grey')
domains <- unique(s2names[s2names$Factor %in% mnames,c('Factor','Color')])
domains$Domain <- color2domain[domains$Color]
rownames(domains) <- domains$Factor
# Compute distances
#source('../R/distGPS.R')
d.s2 <- distGPS(GRangesList(s2.repset),metric='avgdist',mc.cores=12)
d.bg3 <- distGPS(GRangesList(bg3.repset),metric='avgdist',mc.cores=12)
#s2names.repset <- unique(s2names[s2names$Factor %in% mnames,c('Factor','Color','Domain')])
#bg3names.repset <- unique(bg3names[bg3names$Factor %in% mnames,c('Factor','Color','Domain')])
# Compute inter-domain distances
dd.s2 <- domainDist(as.matrix(d.s2),gps='factors',domain=domains$Color,type='inter')
dd.bg3 <- domainDist(as.matrix(d.bg3),gps='factors',domain=domains$Color,type='inter')
# Random seed
set.seed(149)
# Alterate s2
s2.alt <- s2.repset
s2.alt[['EZ']] <- GRanges(rbind(as.data.frame(s2.repset[['GAF']])[sample(1:(length(s2.repset[['GAF']])/2)),],as.data.frame(s2.repset[['HP1B']])[sample(1:(length(s2.repset[['HP1B']])/2)),]))
d.s2.alt <- distGPS(GRangesList(s2.alt),metric='avgdist',mc.cores=12)

# Plot S2 vs BG3
par(las=1,mar=c(4,8,4,4))
mycors1 <- rev(diag(cor(as.matrix(d.s2),as.matrix(d.bg3))))
barplot(mycors1,horiz=TRUE,xlim=c(0,1),main='S2 / BG3',col=domains[names(mycors1),'Color'],font=2)
for (i in 1:length(summary(mycors1))) abline(v=summary(mycors1)[i],col=i,lwd=2,lty=3)

# Plot S2 Altered vs BG3
par(las=1,mar=c(4,8,4,4))
mycors2 <- rev(diag(cor(as.matrix(d.s2.alt),as.matrix(d.bg3))))
barplot(mycors2,horiz=TRUE,xlim=c(0,1),main='S2 Altered / BG3',col=domains[names(mycors2),'Color'],font=2)
for (i in 1:length(summary(mycors2))) abline(v=summary(mycors2)[i],col=i,lwd=2,lty=3)

###################################################
### code chunk number 47: minimum factor set
### https://github.com/singlecoated/chroGPS2/tree/master/examples/rankFactors
###################################################

source('../R/rankFactors.R')
# Call rankFactors for HP1a repression domain, select a combination of 4 factors
library(caTools)
rank.factors.4 <- rankFactors(d.s2,domains,ranktype='domainDist',selName='Domain',selValue='HP1a',k=3,mc.cores=4)
ddd <- as.data.frame(do.call(rbind,lapply(rank.factors.4,unlist)))
ddd <- ddd[order(ddd$intra,decreasing=FALSE),]

head(ddd)

# Now plot superimposed results sorted by intra-domain distance
plot(ddd$intra,type='l',col='black',lwd=6,axes=FALSE,xlab='',ylab='',main='HP1a repression')
axis(2)
par(new=TRUE)
plot(ddd$inter,pch=19,col='grey',axes=FALSE,xlab='',ylab='')
axis(4)
par(new=TRUE)
plot(smooth.spline(ddd$inter),type='l',col='darkgrey',lwd=6,axes=FALSE,xlab='',ylab='')
legend('topright',lwd=6,col=c('black','grey'),legend=c('IntraDomain','InterDomain'),bg='white')
axis(1)

# Perform computation
glm.rank <- rankFactorsGlobally(s2.tab,ranktype='glm',glm.threshold=0.75,mc.cores=12)

# Returned objects are lists named by the factor with highest prediction accuracy in each interation, as well as the rest, let's generate a matrix
library(gtools)
glm.rank <- do.call(smartbind,glm.rank)

# Now let's order based on which factor is removed in each iteration
glm.rank <- glm.rank[,rownames(glm.rank)]
glm.rank[1:5,1:5]

# And finally some plots...
boxplot(glm.rank,horizontal=TRUE,lwd=2,ylab='',xlab='Prediction accuracy')

###################################################
### code chunk number 48: region-based map
### https://github.com/singlecoated/chroGPS2/tree/master/examples/diffPromoters
###################################################

# Define promoters as 1kb upstream of TSS
# + strand, 1kb <<<
# - strand, 1kb >>>
promoters <- annot
end(promoters) <- ifelse(promoters$strand==1,start(promoters),end(promoters)+1000)
start(promoters) <- end(promoters)-1000

# Intersect, DEPRECATED, now we keep things over promoters
s2 <- s2raw
#s2.i <- RangedDataList(lapply(as.list(s2),function(x) do.call(rbind,mclapply(names(x),function(chr) RangedData(intersect(ranges(x)[[chr]],ranges(promoters)[[chr]]),space=chr),mc.cores=8,mc.preschedule=FALSE))))
s2.i <- RangedDataList(mclapply(as.list(s2),function(x) x[x %over% promoters,],mc.cores=12))

# Make map
d <- distGPS(s2.i,metric='avgdist',mc.cores=12)
m <- mds(d,type='isoMDS')
fnames <- s2names$Factor
s2names$Color[s2names$Color=='grey'] <- 'orange'
fcolors <- paste(col2hex(s2names$Color),'BB',sep='')
bcolors <- paste(col2hex(s2names$Color),'FF',sep='')

# Revert Y points to match orientation of original S2
m@points[,2] <- -1 * m@points[,2]

# Plot
plot(m,drawlabels=TRUE,labels=s2names$Factor,point.pch=19,point.cex=4,text.cex=0.75,point.col=s2names$Color,main='S2 @ Promoters (Avgdist-isoMDS)',font=2,xlim=c(-1,1),ylim=c(-1,1))
legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',m@R.square,m@stress))
legend('topright',legend='',bty='n')
\end{lstlisting}

###################################################
### code chunk number 49: diff factors maps
### https://github.com/singlecoated/chroGPS2/tree/master/examples/diffOrigins
###################################################

# Load RepliSeq data (used for Aleix H2av)
orig <- system('ls /Volumes/biostats/consulting/ferran_azorin/jbernues_brdu/reports/tables/files_4_igv/RepliSeq*.bed',intern=TRUE)
names(orig) <- c('Early','Early.Mid','Late','Late.Mid')
orig <- lapply(orig,read.delim,sep='\t',header=FALSE,as.is=TRUE)
orig <- lapply(orig,function(x) { colnames(x) <- c('space','start','end'); RangedData(x) })
orig <- orig[c('Early','Early.Mid','Late.Mid','Late')]

# Intersect s2 with repliseq, filter peaks by overlap with origin set
s2 <- s2raw
s2.origs <- lapply(orig,function(o) RangedDataList(mclapply(as.list(s2),function(x) x[x %over% o,],mc.cores=8)))

# Make distance sets
d.origs <- lapply(s2.origs,function(x) distGPS(x,metric='avgdist',mc.cores=12))
m.origs <- lapply(d.origs,mds,type='isoMDS')

# Modify colors and add some transparency
fnames <- s2names$Factor
s2names$Color[s2names$Color=='grey'] <- 'orange'
fcolors <- paste(col2hex(s2names$Color),'BB',sep='')
bcolors <- paste(col2hex(s2names$Color),'FF',sep='')

m1 <- m.origs[['Early-Mid']]
 m2 <- m.origs[['Late']]

## Perform differential Procrustes analysis
pp <- procrustesDiff(m1,m2)

## Plot both maps before and after adjustment
m3 <- new('mds',points=pp$Yrot)
plot(0,xlim=c(-1,1),ylim=c(-1,1),xlab='',ylab='',xaxt='n',yaxt='n',col='NA')
segments(m1@points[,1],m1@points[,2],m3@points[,1],m3@points[,2],col='red')
par(new=TRUE)
plot(m1,drawlabels=TRUE,labels=s2names$Factor,point.pch=19,point.cex=4,text.cex=0.75,point.col=s2names$Color,main=sprintf('S2@Origins, adjusted (Avgdist-isoMDS)'),font=2,xlim=c(-1,1),ylim=c(-1,1))
par(new=TRUE)
plot(m3,drawlabels=TRUE,labels=s2names$Factor,point.pch=19,point.cex=4,text.cex=0.75,point.col=s2names$Darkcolor,text.col='grey',main='',xaxt='n',yaxt='n',font=2,xlim=c(-1,1),ylim=c(-1,1))

## Plot Procrustes errors
plot(pp)
par(las=1,mar=c(4,12,4,4)); barplot(sort(residuals(pp),decr=TRUE),horiz=TRUE,xlim=c(0,max(residuals(pp))+.1),col=heat.colors(length(residuals(pp))),main='Procrustes errors')


###################################################
### code chunk number 50: diff species maps
### https://github.com/singlecoated/chroGPS2/tree/master/examples/diffSpecies
###################################################

# Compute join and separate distances
load(file=file.path(rdatadir,'K562_clean_distGPS_Avg.RData'))
d1 <- d

## MDS
m.iso <- mds(d,type='isoMDS')
m1 <- m.iso # we choose isoMDS solution
m1@points <- scale(m1@points) / max(scale(m1@points))

# Load sampleinfo
sampleinfo <- read.csv(file.path(datadir,'sampleinfo.csv'),row.names=1,as.is=TRUE)
sampleinfo <- sampleinfo[sampleinfo$filesize!=0,]
sampleinfo$factor <- toupper(unlist(lapply(strsplit(sampleinfo$antibody,'_'),function(x) x[1])))

# Load S2 data, build map and apply procrustes
load('/Volumes/biostats/research/chroGPS/data/Paper/Main/s2_Jan13.RData')
load('/Volumes/biostats/research/chroGPS/reports/figs/Paper/chroGPS_v9/F1.RData')
d2 <- d
fnames <- s2names$Factor
s2names$Color[s2names$Color=='grey'] <- 'orange'
fcolors <- paste(col2hex(s2names$Color),'BB',sep='')
bcolors <- paste(col2hex(s2names$Color),'FF',sep='')
m2 <- mds(d,type='isoMDS')

# First Procrustes, anchor K562, adjustS2, it is the same anyway
# Join mds solutions and distance matrixes
m <- new('mds',points=rbind(m1@points,m2@points),Type='isoMDS',Adj=FALSE,R.square=-1,stress=-1)
d <- matrix(NA,nrow=nrow(d1@d)+nrow(d2@d),ncol=nrow(d1@d)+nrow(d2@d))
colnames(d) <- rownames(d) <- c(rownames(m1@points),rownames(m2@points))
d[1:nrow(d1@d),1:ncol(d1@d)] <- d1@d
d[(nrow(d1@d)+1):(nrow(d1@d)+nrow(d2@d)),(nrow(d1@d)+1):(nrow(d1@d)+nrow(d2@d))] <- d2@d
d <- new('distGPS',d=d,metric='avgdist',type='factors')

# Apply procrustes
id1 <- toupper(sampleinfo$factor)
id2 <- toupper(s2names$Factor)
intersect(id1,id2)
mm <- procrustesAdj(m,d,adjust=c(rep('K562',nrow(d1@d)),rep('S2',nrow(d2@d))),sampleid=c(id1,id2))

# Use S2 colors on K562 factors
colorinfo <- unique(s2names[,c('Factor','Color')]); rownames(colorinfo) <- colorinfo$Factor
sampleinfo$Color <- colorinfo[sampleinfo$factor,'Color']
sampleinfo$Color[grep('POL2',sampleinfo$factor)] <- 'purple'
sampleinfo$Color[grep('EZ',sampleinfo$factor)] <- 'yellow'
sampleinfo$Color[is.na(sampleinfo$Color)] <- 'lightgrey'

setwd(figsdir)

# Plot
pdf('K562_clean_S2_oct17.pdf',height=12,width=12)
# Colors and shapes
mypchs=c(rep(19,nrow(d1@d)),rep(15,nrow(d2@d)))
mycols=c(sampleinfo$Color,s2names$Color)
# Labels
labels <- c(sampleinfo$factor,s2names$Factor)
mytextcols <- ifelse(labels %in% sampleinfo$factor & labels %in% s2names$Factor,'black','grey')
# Human
#plot(m1,drawlabels=TRUE,labels=sampleinfo$factor,point.pch=19,point.cex=4,text.cex=0.75,text.col=mytextcols,point.col=sampleinfo$Color,main='Human K562 (Avgdist-isoMDS)',font=2)
#legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',m1@R.square,m1@stress))
# Human, zoom in
#plot(m1,drawlabels=TRUE,labels=sampleinfo$factor,point.pch=19,point.cex=4,text.cex=0.75,text.col=mytextcols,point.col=sampleinfo$Color,main='Human K562 (Avgdist-isoMDS)',font=2,xl#im=c(-0.5,0.5),ylim=c(-0.5,0.5))
#legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',m1@R.square,m1@stress))
# Fly
#plot(m2,drawlabels=TRUE,labels=s2names$Factor,point.pch=19,point.cex=4,text.cex=0.75,text.col=mytextcols,point.col=s2names$Color,main='Dmel S2 (Avgdist-isoMDS)',font=2)
#legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',m2@R.square,m2@stress))
# Fly, zoom in
#plot(m2,drawlabels=TRUE,labels=s2names$Factor,point.pch=19,point.cex=5,text.cex=0.75,text.col=mytextcols,point.col=s2names$Color,main='Dmel S2 (Avgdist-isoMDS)',font=2,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
#legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',m2@R.square,m2@stress))
# Procrustes
plot(mm,drawlabels=TRUE,labels=labels,point.pch=mypchs,point.cex=5,text.cex=0.9,text.col=mytextcols,point.col=mycols,main='Human K562 / Dmel S2 (Avgdist-isoMDS-Procrustes)',font=2,xlim=range(mm@points),ylim=range(mm@points))
ctcf1 <- grep('CTCF',labels)[grep('CTCF',labels)<77]
ctcf2 <- grep('CTCF',labels)[grep('CTCF',labels)>76]
segments(x0=-.12,y0=.35,x1=0.4,y1=0.20,lwd=5,col='red')
mmm <- mm
sel <- mytextcols=='black'
mmm@points <- mmm@points[sel,]
par(new=TRUE)
plot(mmm,drawlabels=TRUE,labels=labels[sel],point.pch=mypchs[sel],point.cex=5,text.cex=0.9,text.col=mytextcols[sel],point.col=mycols[sel],main='',,xaxt='n',yaxt='n',font=2,xlim=range(mm@points),ylim=range(mm@points))
legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f\nMatch points: %s',mm@R.square,mm@stress,paste(intersect(id1,id2),collapse=',')),bty='n')
legend('bottomleft',legend=c('Human','Fly'),pch=c(19,15),pt.cex=4,col='darkgrey',cex=2)
#plot(mm,drawlabels=TRUE,labels=labels,point.pch=mypchs,point.cex=5,text.cex=0.9,text.col=mytextcols,point.col=mycols,main='Human K562 / Dmel S2 (Avgdist-isoMDS-Procrustes)',font=2,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
#plot(mm,drawlabels=TRUE,labels=labels,point.pch=mypchs,point.cex=5,text.cex=0.9,text.col=mytextcols,[M J3point.col=mycols,main='Human K562 / Dmel S2 (Avgdist-isoMDS-Procrustes)',font=2,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
dev.off()

###################################################
### code chunk number 51: diff genes
### https://github.com/singlecoated/chroGPS2/tree/master/examples/diffGenes
###################################################

# Load data
load(file.path(datadir,'s2_Jan13.RData'))
load(file.path(datadir,'bg3_Jan13.RData'))

# Load new function
source('/Volumes/biostats/routines/R/chroGPS2_library/functions/diffGenes.R')

# Summarize factor replicates with method 'any' so that 1 replicate having the mark is enough
s2.tab <- mergeReplicates(s2.tab,s2names$Factor,'any')
bg3.tab <- mergeReplicates(bg3.tab,bg3names$Factor,'any')

# Join, use common factors. Then use common genes only from those that have at least one mark in both s2 and bg3
x <- diffGPS(s2.tab,bg3.tab,'S2','BG3')

# Build map and cluster as always
d <- distGPS(x,metric='tanimoto',uniqueRows=TRUE)
m <- mds(d,type='classic',splitMDS=TRUE,split=0.16,mc.cores=4)
mm <- mds(d,m,type='boostMDS',samplesize=0.005,mc.cores=6) # Do this in pipeline or Aurora please, it stresses out memory quite a bit. review a final warning
unique(rownames(m@points)==rownames(mm@points)) # sanity check
# This should be incorporated in the function code...
m@points <- m@points[rownames(d@d),]
mm@points <- mm@points[rownames(d@d),]
save(m,mm,file=file.path(diffdatadir,'mds_S2_BG3_tanimoto_classic_split_Jan18.RData'))

# Cluster
h <- hclust(as.dist(d@d),method='average')
#source('/Volumes/biostats/routines/R/chroGPS_library/chroGPS/R/clusGPS.R')
clus <- clusGPS(d,mm,h,k=max(cutree(h,h=0.5)),ngrid=10000,mc.cores=8,recalcDist=FALSE,verbose=FALSE)
clus.merged <- mergeClusters(clus,brake=0,mc.cores=8)
save(clus,clus.merged,file=file.path(diffdatadir,'clusGPS_merged_S2_BG3_tanimoto_weighted_classic_boost_hclust_avg_Jan18.RData'))

# Use new function to profile clusters
source('/Volumes/biostats/routines/R/chroGPS2_library/functions/profileClusters.R')
pc <- profileClusters2(x,clus.merged,normalize=TRUE)
#pc <- profileClusters2(x,clus,normalize=TRUE) # non/merged clusters
heatmap.2(pc,trace='none',scale='none',col=bluered(100))

# Perform 'differential' analysis
x.diff <- res <- diffGPS.clus(x,mm,clus.merged,label.x='S2',label.y='BG3',clusName=clusNames(clus.merged)[1],fdr=TRUE,mc.cores=8)
write.csv(x.diff,'kk_fdrest.csv')

# Select genes changing clusters with FDR 0.05
xx.diff <- x.diff[x.diff$ClusID.S2!=x.diff$ClusID.BG3 & x.diff$FDR.S2<0.25 & x.diff$FDR.BG3<0.25,]
xx.diff$CC <- paste(xx.diff$ClusID.S2,xx.diff$ClusID.BG3,sep='.')
table(xx.diff$CC)
write.csv(xx.diff,'kk_fdrest2.csv')

# Perform enrichment test using getEnrichedGO from chippeakanno package
library(ChIPpeakAnno)
enriched.GO <- lapply(sort(unique(xx.diff$CC)),function(cc) {
        fbid <- as.character(xx.diff$geneid[xx.diff$CC==cc])
            if (length(fbid)>=25)
                        ans <- getEnrichedGO(annotatedPeak=fbid,orgAnn='org.Dm.eg.db',maxP=0.05,multiAdj=TRUE,multiAdjMethod='BH')
            else ans <- NULL
            return(ans)
    })
names(enriched.GO) <- sort(unique(xx.diff$CC))
enriched.GO <- enriched.GO[unlist(lapply(enriched.GO,length))>0]
enriched.GO <- lapply(enriched.GO,function(x) lapply(x,function(y) unique(y[,-ncol(y)])))

# Plot results with diffGPS.plot function
source('/Volumes/biostats/routines/R/chroGPS2_library/functions/diffFactors.R')
source('/Volumes/biostats/routines/R/chroGPS2_library/functions/diffGenes.R')
res.sel <- res[res$ClusID.S2!=res$ClusID.BG3,]

diffGenes.plot(x,mm,clus.merged,res.sel,transitions='10.2',label.x='S2',label.y='BG3',fdr1=0.25,fdr2=0.25)

# Now plot enriched GO info
txt <- paste("allsites <- c(",paste("s2[[",c(33,44),"]]",collapse=','),")",sep='')
eval(parse(text=txt))

# Test mergeReplicates for intervals directly
source('/Volumes/biostats/routines/R/chroGPS2_library/functions/mergeReplicates.R')
kk1 <- mergeReplicateList(s2,id=s2names$Factor,mergeBy='any',mc.cores=8)
kk2 <- mergeReplicateList(s2,id=s2names$Factor,mergeBy='all',mc.cores=8)
kk3 <- mergeReplicates(s2,id=s2names$Factor,mergeBy=0.3,mc.cores=8)
kk4 <- mergeReplicates(s2,id=s2names$Factor,mergeBy=1,mc.cores=8)
mmm <- mergeReplicates(m.factors,id=s2names$Factor)

# Permutation tests...

