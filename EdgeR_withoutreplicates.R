
bioclimEdgeR<- function(data,biodata, bioclim, filtercpm){
library(edgeR) # version 3.16.5
#load counts
#"/home/fbalao/Datos/ARTICULO/Dacthylorhizapolyploids/Dactylorhiza_4xExpression/data/19.01.2017.counts.txt"
counttable = read.table(data, header = TRUE, row.names = 1) 
colnames(counttable)[1:4]<-c("tB1830_1", "tB1830_2", "tB1833_1", "tB1833_2") #Pool ebudensis and traunsteineri
head(counttable)

# #Load climatic variables
#load(biodata) #This load the bioclim variables
#rownames(biodata2)<-colnames(counttable)
biodata3 <- as.data.frame(biodata)


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
 
 #filter (optional?)
 keep <- rowSums(cpm(e)>filtercpm) >= 9
 e <- e[keep, , keep.lib.sizes=FALSE]
 
 
 #recompute library sizes 
 e$samples$lib.size <- colSums(e$counts)

#Normalizing
e <- calcNormFactors(e,method = "upperquartile")
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

#' Plot Gene Expression vs Enviromental data
#'
#' @param bioclim Enviromental variable (bioclim_1:20)
#' @param cluster Contig to plot
#' @param counttable Count data
#'
#' @return
#' @export
#'
#' @examples plotbio("bio_11","Cluster-38451.7", counttable3)
plotbio<-function(bioclim,cluster,normcounttable){
colspecies<-factor(substr(colnames(normcounttable),1,1))
plot(biodata4[bioclim][,1], cpm(normcounttable$counts)[cluster,] , cex=1, pch=16, ylab="Normalized counts (CPM)", xlab=bioclim, col=c(2,3)[colspecies], main="")
legend("topright",legend = c("D. majalis", "D. traunsteineri"), col=c(2,3), pch=16, cex=0.7)
}

####### Loop
# res<-list()
# var<-colnames(biodata2)
# for (i in 1: length(colnames(biodata2))){
#   res[[i]]<-bioclimEdgeR(var[i])
#   }
# names(res)<-var