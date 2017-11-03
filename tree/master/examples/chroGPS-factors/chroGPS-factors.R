> load('/Volumes/biostats/research/chroGPS/data/Paper/Main/s2_Jan13.RData')

> class(s2)
[1] "RangedDataList"
attr(,"package")
[1] "IRanges"

> head(names(s2))
[1] "ASH1-Q4177.S2"    "BEAF-70.S2"       "BEAF-HB.S2"       "Chro(Chriz)BR.S2"
[5] "Chro(Chriz)WR.S2" "CP190-HB.S2"     

> s2
RangedDataList of length 76
names(76): ASH1-Q4177.S2 BEAF-70.S2 BEAF-HB.S2 ... WDS_Q2691.S2 ZW5.S2

> d <- distGPS(s2,metric='avgdist',mc.cores=8)
>
> d
Object of class distGPS with avgdist distances between 76 objects 
>

> as.matrix(d)[1:5,1:4]
                 ASH1-Q4177.S2 BEAF-70.S2 BEAF-HB.S2 Chro(Chriz)BR.S2
ASH1-Q4177.S2        0.0000000  0.9302190  0.7067805       0.71324407
BEAF-70.S2           0.9302190  0.0000000  0.3450349       0.34739217
BEAF-HB.S2           0.7067805  0.3450349  0.0000000       0.13782518
Chro(Chriz)BR.S2     0.7132441  0.3473922  0.1378252       0.00000000
Chro(Chriz)WR.S2     0.7320215  0.3479041  0.1376213       0.04270634
> 

> m <- mds(d,k=2,type='classic')

> m
Object of class MDS approximating distances between 76 objects 
R-squared= 0.7402 Stress= 0.1678 

> m <- mds(d,k=3,type='classic')

> m
Object of class MDS approximating distances between 76 objects 
R-squared= 0.8275 Stress= 0.0888 

> m.iso <- mds(d,k=2,type='isoMDS')

> m.iso
Object of class MDS approximating distances between 76 objects 
R-squared= 0.8294 Stress= 0.0737 
> 
>
> head(getPoints(m.iso))
                        [,1]        [,2]
ASH1-Q4177.S2    -0.42154547 -0.32763008
BEAF-70.S2       -0.42089405  0.26956091
BEAF-HB.S2       -0.18283057  0.09523077
Chro(Chriz)BR.S2 -0.09120635  0.07736748
Chro(Chriz)WR.S2 -0.10749822  0.09335583
CP190-HB.S2       0.04427325 -0.04324142
>
> head(s2names)[,-4]
    ExperimentName      Factor             ID      Color Darkcolor
1    ASH1-Q4177.S2        ASH1 modENCODE_2984 lightgreen darkgreen
2       BEAF-70.S2        BEAF  modENCODE_922       grey  darkgrey
3       BEAF-HB.S2        BEAF  modENCODE_274       grey  darkgrey
4 Chro(Chriz)BR.S2 CHRO(CHRIZ)  modENCODE_278 lightgreen darkgreen
5 Chro(Chriz)WR.S2 CHRO(CHRIZ)  modENCODE_279 lightgreen darkgreen
6      CP190-HB.S2       CP190  modENCODE_925       grey  darkgrey
>
plot(m.iso,point.col=s2names$Color,point.cex=8,text.cex=.75,drawlabels=TRUE,labels=s2names$Factor,font=2)
legend('topleft',legend=sprintf('R2=%.3f / Stress=%.3f',getR2(m.iso),getStress(m.iso)))
