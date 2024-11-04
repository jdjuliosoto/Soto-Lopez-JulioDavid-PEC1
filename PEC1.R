# Conectarse a github
library(usethis)
usethis::create_github_token()
library(gitcreds)
gitcreds::gitcreds_set()

# Data set: Gastric_cancer

# Cargando paquetes
library(readxl)
library(Biobase)
library(SummarizedExperiment)
library(usethis)
library(gitcreds)

# recuperando los datos 

GastricCancer_NMR <- read_excel("GastricCancer_NMR.xlsx")
Peak <- read_excel("GastricCancer_NMR.xlsx", sheet = "Peak")

Peak$filtrado <- subset(Peak[,4] <= 10 & Peak[,5] <= 20) #filtrando
Peak_filtrado <- Peak[Peak$filtrado == TRUE,]
GastricCancer <- GastricCancer_NMR[GastricCancer_NMR$Idx %in% 
                                     Peak_filtrado$Idx,]

# DE .xlsx A SummarizedExperiment

# datos de la expresion a matriz traspuesta
Expresiongenes <- t(as.matrix(GastricCancer[,5:153]))
colnames(Expresiongenes) <- paste0("sample_", 
                                   1:48)# agrego nombre a las columnas

# phenodata
targets <- data.frame(sampleNames = paste0("sample_", 1:48),
                      sampleType = GastricCancer$SampleType,
                      sampleIdx = GastricCancer$Idx,
                      sampleClass = GastricCancer$Class,
                      row.names = 1)

columnDesc <- data.frame(labelDescription= c("Sample/QC",
                                             "Identificador",
                                             "Subclass (QC: Control de calidad,
                                             GC: Cancer gastrico,
                                             BN: Tumor benigno,
                                             HE: Control saludable)"))

myAnnotDF <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)

# agrego nombre de genes
mymetabolitos <- paste0("M", 1:149)

# informacion minima del experimento
myInfo=list(myName="Julio Soto",
            myLab="Analisisd e datos omicos",
            myContact="jdjulio@uoc.es",
            myTitle="PEC1")

myDesc <- new("MIAME", name = myInfo[["myName"]],
              lab = myInfo[["myLab"]],
              contact = myInfo[["myContact"]] ,
              title = myInfo[["myTitle"]],
              url = "https://cimcb.github.io/MetabWorkflowTutorial/Tutorial1.html",
              abstract = "Dataset used in the CIMBC tutorial on Basic 
              Metabolomics Data Analysis Workflow 
              The tutorial describes the data as follows:
              - The study used in this tutorial has been previously published 
              as an open access article Chan et al. (2016), in the British 
              Journal of Cancer. - The deconvolved and annotated data file have 
              been deposited at the Metabolomics Workbench data repository 
              (Project ID PR000699). - The data can be accessed directly via 
              its project DOI:10.21228/M8B10B - 1H-NMR spectra were acquired 
              at Canada’s National High Field Nuclear Magnetic Resonance Centre 
              (NANUC) using a 600 MHz Varian Inova spectrometer. - Spectral 
              deconvolution and metabolite annotation was performed using 
              the Chenomx NMR Suite v7.6. ")

# Datos a tipo ExpressionSet
myEset <- ExpressionSet(assayData = Expresiongenes,
                        phenoData = myAnnotDF,
                        featureNames = mymetabolitos,
                        experimentData = myDesc)

# SummarizedExperiment
eset <- makeSummarizedExperimentFromExpressionSet(from = myEset, 
                                                  mapFun = naiveRangeMapper)
# resumen del objeto sumarizedExperiment
class(eset)
show(eset)

# explorando el objeto sumarizedExperiment
assays(eset)$exprs[1:7,1:7]
colData(eset)
metadata(eset)

# estadistica univariada
summary(assays(eset)$exprs[1:7,1:7])

# visualizacion 

colores <- c("red", "blue", "yellow", "green")  
gctotal <- colores[as.factor(targets$sampleClass)]

boxplot(log10(assays(eset)$exprs), 
        main="Concentracion de metabolitos para todas las muestras",
        xlab="Muestra", col = gctotal,
        ylab="Log 10 Concentracion", las=2, cex.axis=0.5, cex.main=0.9)


# visualizacion multivariante

# cambiar los NA por 0
library(tidyverse)
library(VIM)
library(scales)

x <- assays(eset)$exprs
x_ <- log10(x)
x_1 <- assays(eset)$exprs #backup
Xscale <- scale(x_) # escalar

Xknn <- kNN(Xscale, k = 3, imp_var = FALSE) # imputar valores NA
head(Xknn)

apply(Xknn, MARGIN = -1, function(x) sum(is.na(x))) # evaluar na


# PCA
library(FactoMineR)
library(factoextra)
pca <- PCA(X = Xknn, scale.unit = TRUE, graph = FALSE)
head(pca$eig)
fviz_pca_ind(pca, geom.ind = "point", 
             col.ind = "#FC4E07", 
             axes = c(1, 2), 
             pointsize = 1.5) 

# Cluster
clust.euclid.average <- hclust(dist(t(Xknn)),method="average")
plot(clust.euclid.average, hang=-1)
