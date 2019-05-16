###########################################################################################
#               Msc. PELAYO G. DE LENA RODRÍGUEZ                                          #
#               R version 3.5.3 (2019-03-11) 
#                 RStudio Version 1.2.1335
#         Tutorials for the R/Bioconductor Package copynumber
#     Instituto de Estudios Celulares y Moleculares - ICM                                #
#                     MADRID, 16/05/2019                                                 #
###########################################################################################


# start with a clean slate

rm(list=ls(all=TRUE)) 

# INSTALL PACKAGE FROM bIOcONDUCTOR

source("https://bioconductor.org/biocLite.R")
biocLite("copynumber")

# load library

library(copynumber)

# load, inspect and subset dataset
data(lymphoma)
View(lymphoma)
sub.lymphoma <- subsetData(data=lymphoma)
sub.lymphoma[1:10,]


# use winsorize function to detect outliers
lymph.wins <- winsorize(data=sub.lymphoma,verbose=FALSE)
lymph.wins[1:10,]

# retrieve results
wins.res <- winsorize(data=sub.lymphoma,
                      return.outliers=TRUE,
                      verbose=FALSE)

wins.res$wins.outliers[1:10,]

# segmentation algorithm pcf
single.seg <- pcf(data=lymph.wins,gamma=12,verbose=FALSE)

plotGenome(data=sub.lymphoma,
           segments=single.seg,
           sample=1,
           cex=3)


plotSample(data=sub.lymphoma,
           segments=single.seg,
           layout=c(5,5),
           sample=1,
           cex=3)

multi.seg <- multipcf(data=lymph.wins,verbose=FALSE)
head(multi.seg)

# not run
plotChrom(data=lymph.wins,
          segments=multi.seg,
          layout=c(2,1))



# load, inspect the dataset
data(logR)
data(BAF)

# use winsorize function to detect outliers
logR.wins <- winsorize(logR,verbose=FALSE)


# segmentation algorithm
allele.seg <- aspcf(logR.wins,
                    BAF,
                    verbose=FALSE)
head(allele.seg)

plotAllele(logR,
           BAF,
           allele.seg,
           sample=1,
           chrom=c(1:4),
           layout=c(2,2))

lymphoma.res <- pcf(data=lymphoma,
                    gamma=12,
                    verbose=FALSE)

chr.from <- c(2,12,4)
pos.from <- c(168754669,847879349,121809306)
chr.to <- c(14,21,17)
pos.to <- c(6147539,301955563,12364465)
cl <- c(1,1,2)
arcs <- cbind(chr.from, pos.from, chr.to, pos.to,cl)

plotCircle(segments=lymphoma.res,
           thres.gain=0.2,
           arcs=arcs)

plotHeatmap(segments=lymphoma.res,upper.lim=0.3)
plotAberration(segments=lymphoma.res,thres.gain=0.2)







# MODEL PARAMETERS: gamma
data(micma)
View(micma)
plotGamma(micma,chrom=17,cex=3)


