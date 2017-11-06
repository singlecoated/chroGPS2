# Call rankFactors for HP1a repression domain, select a combination of 4 factors
rank.factors.4 <- rankFactorsInDomain(d,domains,selName='Domain',selValue='HP1a repression',k=4,mc.cores=12)
ddd <- as.data.frame(do.call(rbind,lapply(rank.factors.4,unlist)))
ddd <- ddd[order(ddd$intra,decreasing=FALSE),]

> head(ddd)
                                                                                    inter     intra
H3K9me2-Ab2-(new-lot).S2 + H3K9me2-antibody2.S2 + H3K9me3.S2 + Su(var)3-9.S2    0.9241429 0.4058082
H3K9me2-antibody2.S2 + H3K9me3.S2 + HP1a_wa184.S2 + Su(var)3-9.S2               0.9036224 0.4060373
H3K9me2-Ab2-(new-lot).S2 + H3K9me2-antibody2.S2 + H3K9me3.S2 + HP1a_wa184.S2    0.8982983 0.4061356
H3K9me2-Ab2-(new-lot).S2 + H3K9me3.S2 + HP1a_wa184.S2 + Su(var)3-9.S2           0.9135723 0.4062027
H3K9me2-antibody2.S2 + H3K9me3.S2 + HP1a_wa191.S2 + Su(var)3-9.S2               0.8984177 0.4062145
H3K9me2-Ab2-(new-lot).S2 + H3K9me2-antibody2.S2 + HP1a_wa184.S2 + Su(var)3-9.S2 0.9126298 0.4062168

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

> # Basic data
> head(s2.tab[,1:5])
            ASH1-Q4177.S2 BEAF-70.S2 BEAF-HB.S2 Chro(Chriz)BR.S2 Chro(Chriz)WR.S2
FBgn0000003             0          0          0                0                0
FBgn0000008             0          0          0                1                1
FBgn0000014             0          0          0                0                0
FBgn0000015             0          0          0                0                0
FBgn0000017             1          0          1                1                1
FBgn0000018             0          0          1                1                1

> # Perform computation
> glm.rank <- rankFactorsGlobally(s2.tab,ranktype='glm',glm.threshold=0.75,mc.cores=12)

> # Returned objects are lists named by the factor with highest prediction accuracy in each interation, as well as the rest, let's generate a matrix
> glm.rank <- do.call(smartbind,glm.rank)

# Now let's order based on which factor is removed in each iteration
> glm.rank <- glm.rank[,rownames(glm.rank)]
> glm.rank[1:5,1:5]
                    EZ.Q3421.S2 Su.var.3.9.Q2598.S2     Ez.S2 Su.Hw..VC.S2 HP2..Ab2.90..S2
EZ.Q3421.S2           0.9960983           0.9800479 0.9796045    0.9768555       0.9760575
Su.var.3.9.Q2598.S2          NA           0.9800479 0.9789838    0.9765008       0.9760575
Ez.S2                        NA                  NA 0.9790724    0.9765008       0.9757028
Su.Hw..VC.S2                 NA                  NA        NA    0.9765008       0.9756141
HP2..Ab2.90..S2              NA                  NA        NA           NA       0.9756141

# And finally some plots...
> boxplot(glm.rank,horizontal=TRUE,lwd=2,ylab='',xlab='Prediction accuracy')
