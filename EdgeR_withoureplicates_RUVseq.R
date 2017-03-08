bioclimEdgeR_RUVseq<- function(bioclim){
  library(edgeR) # version 3.16.5
#load counts
data = read.table("/home/fbalao/Datos/ARTICULO/Dacthylorhizapolyploids/Dactylorhiza_4xExpression/data/19.01.2017.counts.txt", header = TRUE, row.names = 1) 
source("./RUVseqModel.R") #load RUV stuff
load("~/Datos/ARTICULO/Dacthylorhizapolyploids/ResultsRUV_bio3.RData") #This load the bioclim variables
colnames(data)[1:4]<-c("tB1830_1", "tB1830_2", "tB1833_1", "tB1833_2") #Pool ebudensis and traunsteineri


dat = selectSpecies(data, "t", "m", " ", " ", " ")

rownames(biodata2)<-colnames(data)
biodata3 = biodata2[colnames(dat),]
biodata3 <- as.data.frame(biodata3)
# you can also choose "f" and "i"for fuchsii or incarnata if wanted in analysis.
# furthermore, if you are interested only in pyrenees then "P", alps "A", scan "S" or britain "B"
# eg. selectSpecies(data, "f", "m", "A", "P", " ") means
# all fuchsii and all majalis plus all the pyrenees samples and the alps samples, useful for special comparisons.

set = makeRUVset(dat = dat)
# formats the data so that it can be normalised

normSet = makeRUVrepNormalisedSet(set = set, empircalNonDEgenes = rownames(counts(set)), numFactors = 2, replicates = batchMatrix)
# normalises the data by removing numFactors of unwanted variance. batchMatrix was created in selectSpecies
# normSet is and S4 object, it contains the normalizedCounts (unfortunately not directly usable) & 
# W (W_1 ... W_numFactors) the factors of unwanted variance


# This is fitting a non-linear regression e.g cuadratic, cubic etc
# require(splines)
# envvar<-biodata3$bio_5# Change the climatic variable
# X <- ns(envvar, df=3) # Any df between 3 - 5 usually works well.
#popfactor<-factor(substr(rownames(biodata3),3,6)) #make population factor

markreplicates<-c(2,4,10,14,19,21,23,27)
normSet2<-normSet
normSet2$W<-normSet2$W[-markreplicates,]
normSet2$normalizedCounts<-normSet2$normalizedCounts[,-markreplicates]
set2<-newSeqExpressionSet(counts=counts(set)[,-markreplicates])
env<-as.vector(biodata3[bioclim][-markreplicates,1])

findEdgeRgenes(model = ~env+normSet2$W, group = NULL, set = set2, pvalue = 0.05)
# As you see, the model for DE takes normSet$W correction factors into account. Also, set = set, edgeR never uses normalizedCounts.
# If you change biodata3$bio_3 by X we run a non-linear model
# We can also include a 
hist(resEdgeRtopTags$table$PValue, breaks = 100)
hist(resEdgeRtopTags$table$FDR, breaks = 100)
hist(resEdgeRtopTags$table$logFC, breaks = 100)
# some tests I was using to get the normalisation ok

head(resEdgeRtopTags$table, n = 20L)
print(resEdgeRglm)
}
# I put result into these two "global" varianble when doing findEdgeRgenes. Not elegant but useful, might have to change this when i publish the code.

#Visualization
colspecies<-factor(substr(colnames(dat),1,1))
plot(biodata3$bio_5, counts(set)['Cluster-6494.1',] , cex=1, pch=16, ylab="Normalized counts", xlab="bio5", col=c(2,3)[colspecies], main="Cluster-28751.0-bio5")
legend(200,180,legend = c("D. majalis", "D. traunsteineri"), col=c(2,3), pch=16, cex=0.7)

