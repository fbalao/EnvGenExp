
bioclimEdgeR<- function(bioclim){
library(edgeR) # version 3.16.5
#load counts
counttable = read.table("/home/fbalao/Datos/ARTICULO/Dacthylorhizapolyploids/Dactylorhiza_4xExpression/data/19.01.2017.counts.txt", header = TRUE, row.names = 1) 
colnames(counttable)[1:4]<-c("tB1830_1", "tB1830_2", "tB1833_1", "tB1833_2") #Pool ebudensis and traunsteineri
head(counttable)

#Load climatic variables
load("~/Datos/ARTICULO/Dacthylorhizapolyploids/ResultsRUV_bio3.RData") #This load the bioclim variables
rownames(biodata2)<-colnames(counttable)
biodata3 <- as.data.frame(biodata2)


#Remove species and replicates
species<-factor(substr(colnames(counttable),1,1))
toremove<-species=="i"| species=="f"
counttable2<-counttable[,-c(5:14)]
markreplicates<-c(2,4,10,14,19,21,23,27)
counttable3<-counttable2[,-markreplicates] #remove replicates
biodata4<-biodata3[colnames(counttable3),] #select bioclim dataset without replicates and f/i

#EdgeR analysis
formula<-as.vector(biodata4[bioclim])[,1]
 edesignbio <- model.matrix(~formula)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=counttable3)
e <- calcNormFactors(e)
e <- estimateGLMCommonDisp(e, edesignbio)
e <- estimateGLMTrendedDisp(e, edesignbio)
e <- estimateGLMTagwiseDisp(e, edesignbio)

## MDS Plot
plotMDS(e, main="edgeR MDS Plot")

## Biological coefficient of variation plot
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")

## Fit the model, testing the coefficient for the trausteineri vs majalis comparison
efit <- glmFit(e, edesignbio)
efit <- glmLRT(efit, coef=2)

## Make a table of results
etable <- topTags(efit, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]
print(etable)
}


#Visualization

plotbio<-function(bioclim,cluster){
colspecies<-factor(substr(colnames(counttable3),1,1))
plot(biodata4[bioclim][,1], counttable3[cluster,] , cex=1, pch=16, ylab="Normalized counts", xlab="bio5", col=c(2,3)[colspecies], main="")
legend("topright",legend = c("D. majalis", "D. traunsteineri"), col=c(2,3), pch=16, cex=0.7)
}
