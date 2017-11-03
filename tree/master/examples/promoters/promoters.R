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
