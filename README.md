We present chroGPS2, a computational framework for differential analysis of epigenomes. Methods are provided for efficient integration and comparison of data from different condi- tions or biological backgrounds, accounting and adjusting for systematic biases in order to provide an efficient and statistically robust base for differential analysis. We also include funcionalities for general data assessment and quality control prior to comparing maps, such as functions to study chromatin domain conservation between epigenomic backgrounds, to detect gross technical outliers and also to help in the selection of candidate marks for de-novo epigenome mapping.

- Availability: https://github.com/singlecoated/chroGPS2 
- Contact: oscar.reina@irbbarcelona.org

In recent years, we have assisted to an unprecedented increase in availability of epigenomics data related to whole-genome distribution of transcription factors, histone modifications and other DNA binding proteins, which helped to establish a deeper knowledge about chro- matin states and topologically associated domains, and to offer a broad scope of genomic regulation from both a functional and structural point of view. However, there remains a need for efficient tools for visualization, functional analysis and comparison of epigenomics data.

Previously, we developed chroGPS, an R package based on dimensionality reduction techniques (namely Multidimensional Scaling, MDS) to measure and visualize genome-wide associations between epigenetic factors and their relationship to functional genetic elements in low dimensional maps. Now we extend this software and introduce novel features to perform differential analysis of epigenome maps using Procrustes, hierarchical clustering and Bayesian density estimation methods. Additionally, we provide functions for general data assessment, quality control and to help in the selection of epigenetic factors for de-novo epigenome mapping. ChroGPS2 is integrated in Bioconductor, an open-source collection of R packages for computational analysis of omics data. We illustrate our approach using two publicly available datasets containing extensively mapped human and fruit fly epigenomes (https://www.ncbi.nlm.nih.gov/geo/info/ENCODE.html). 

References

chroGPS package: https://www.bioconductor.org/packages/release/bioc/html/chroGPS.html

chroGPS publication: https://academic.oup.com/nar/article/42/4/2126/2437048
