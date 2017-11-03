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
