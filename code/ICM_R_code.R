###########################################################################################
#               Msc. PELAYO G. DE LENA RODRÍGUEZ                                          #
#               R version 3.5.3 (2019-03-11) 
#                 RStudio Version 1.2.1335
#         Tutorials for the R/Bioconductor Package SNPRelate
#     Instituto de Estudios Celulares y Moleculares - ICM                                #
#                     MADRID, 16/05/2019                                                 #
###########################################################################################

'''
Installation of the package SNPRelate
Preparing Data
Data formats used in SNPRelate
Create a GWAS SNP GDS File
Format conversion from PLINK text/binary files
Format conversion from VCF files
Data Analysis
LD-based SNP pruning
Principal Component Analysis (PCA)

'''

# Borra el entorno y las variables. Inicio limpio de sesión.

rm(list=ls(all=TRUE)) 

# Configura el directorio de trabajo. 
# Ha de ser el lugar donde esten los ficheros .bed, .fam y .bim

setwd("C:/Users/Usuario/Desktop")

# Instala las librerias y paquetes requeridos.

# Install the package from Bioconductor repository:
  
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")
BiocManager::install("GWASTools")
BiocManager::install("SeqArray")

install.packages("lmtest")
install.packages("MASS")

library(SeqArray)
library(GWASTools)
library(SNPRelate)
library(gdsfmt)
library(MASS)
## SNPRelate -- supported b

# Lectura ficheros plink

bed.fn <- "C:/Users/Usuario/Desktop/Lugo.bed"
fam.fn <- "C:/Users/Usuario/Desktop/Lugo.fam"
bim.fn <- "C:/Users/Usuario/Desktop/Lugo.bim"
vcf.fn <- "C:/Users/Usuario/Desktop/LugoRecode.vcf"


# iniciar variable on the fly

gdsfile <- "snps.gds"



#Format conversion from PLINK text/binary files

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gdsfile, family=TRUE,
                cvt.chr="int", cvt.snpid="int", verbose=TRUE)


#Format conversion from VCF files

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

#Summary

snpgdsSummary("snps.gds")


'''

Data Analysis
gdsfmt and SNPRelate (high-performance computing R packages for multi-core
symmetric multiprocessing computer architectures) 
to accelerate two key 
computations in GWAS: principal component analysis (PCA)
and relatedness analysis using identity-by-descent (IBD) measures.

'''

# Open the GDS file
genofile <- snpgdsOpen("snps.gds")

#Get population information
pop_code <- scan("tab.txt", what=character())

# set seed for reproducibility purposes
set.seed(1000)

'''
LD-based SNP pruning
It is suggested to use a pruned set of SNPs 
which are in approximate linkage equilibrium with 
each other to avoid the strong influence of SNP clusters in 
principal component analysis and relatedness analysis.

'''

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

names(snpset) # all the SNP selected (~51000SNPS)

head(snpset$chr1)  # snp.id per chr


#Get all selected snp id

snpset.id <- unlist(unname(snpset))

head(snpset.id)



'''
Principal Component Analysis (PCA)
The functions in SNPRelate for PCA include 
calculating the genetic covariance matrix from
genotypes, computing the correlation coefficients
between sample loadings and genotypes for each SNP,
calculating SNP eigenvectors (loadings), 
and estimating the sample loadings of a new dataset
from specified SNP eigenvectors.

'''


# Run PCA

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

tab <- data.frame(sample.id = pca$sample.id,
EV1 = pca$eigenvect[,1],    # the first eigenvector
EV2 = pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)

head(tab)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


#Plot the principal component pairs for the first
#four PCs:
dev.off()  
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:5], col=tab$pop, labels=lbls)


#Parallel coordinates plot for the top principal components:
  
dev.off()
datpop <- factor(pop_code)[match(pca$sample.id, tab$sample.id)]
parcoord(pca$eigenvect[,1:30], col=datpop)

#calculate the SNP correlations (CORR) between eigenvectors and SNP genotypes:
  
# Get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)






