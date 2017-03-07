#############
##LIBRARIES##
#############

library("RUVSeq")
library("edgeR")
library("RColorBrewer")
library("cgwtools")

#red <- transp("red",.4)
red = rgb(255, 0, 0, max = 255, alpha = 125, names = "red")
black <- rgb(0, 0, 0, max = 255, alpha = 125, names = "black")
## GLOBAL VARIABLES & FUNCTION FOR DATA ORGANISATION
data = read.table("/home/fbalao/Datos/ARTICULO/Dacthylorhizapolyploids/Dactylorhiza_4xExpression/data/19.01.2017.counts.txt", header = TRUE, row.names = 1)
## Combine replicate

combineRep = function(data){
  data = transform(data,  tB1830=tB1830_1+tB1830_2, tB1833=tB1833_1+tB1833_2, ## ebudensis
                   fP1707=fP1707_1+fP1707_2, ##fuchsii
                   mA1567=mA1567_1+mA1567_2,  mP1722= mP1722_1+mP1722_2, mP1748=mP1748_1+mP1748_2, mS1765=mS1765_1+mS1765_2, ##majalis
                   tB1798=tB1798_1+tB1798_2, tB1805=tB1805_1+tB1805_2, tB1812=tB1812_1+tB1812_2, tS1920=tS1920_1+tS1920_2) ##traunteineri
  data = data[, -grep("_", colnames(data))]
  return(data)
}


# ----+ functions: makeVectors
makeSpeciesVector = function(data){
  indm = with(data, grepl("m", colnames(data))); indm[indm==1] = 2
  indt = with(data, grepl("t", colnames(data))); indt[indt==1] = 3
  inde = with(data, grepl("e", colnames(data))); inde[inde==1] = 3
  indf = with(data, grepl("f", colnames(data))); indf[indf==1] = 5
  indi = with(data, grepl("i", colnames(data))); indi[indi==1] = 6 #modifie here if you want to make ebu => traun
  ind = indm+indt+inde+indf+indi; ind = sub(2, "majalis", ind); ind = sub(3, "traunsteineri", ind); ind = sub(5, "fuchsii", ind); ind = sub(6, "incarnata", ind)
  species <<- ind
  return(ind)
}
makeGeographiesVector = function(data){
  indA = with(data, grepl("A", colnames(data))); indA[indA==1] = 2
  indP = with(data, grepl("P", colnames(data))); indP[indP==1] = 3
  indB = with(data, grepl("B", colnames(data))); indB[indB==1] = 4
  indS = with(data, grepl("S", colnames(data))); indS[indS==1] = 5
  ind = indA+indP+indB+indS; ind = sub(2, "alps", ind); ind = sub(3, "pyrenees", ind); ind = sub(4, "britain", ind); ind = sub(5, "scandinavia", ind)
  #geographies <<- ind
  return(ind)
}

makeBatchMatrix = function(data){
  batchIndexVector = with(data, grepl("_", colnames(data)))
  batchVec = which(batchIndexVector == TRUE)
  return(matrix(batchVec, byrow = TRUE, ncol = 2))
}


#########
## MAPS##
#########
geoMap = c("alps"=17, "britain" = 19, "pyrenees" = 15, "scandinavia" = 18)
colorMap = c("ebudensis"="darkgreen","fuchsii"="purple","incarnata"="orange","majalis"="red","traunsteineri"="blue")


#####################
## HELPER FUNCTIONS##
#####################

# ----+ functions: selectFac
selectSpeciesFac = function(ind){
  return(as.factor(species[ind]))
}
selectGeographiesFac = function(ind){
  return(as.factor(geographies[ind]))
}
selectOriginsFac = function(ind){
  return(as.factor(origins[ind]))
}

##############
## FUNCTIONS##
##############


# ----+ function:  selectSpecies
selectSpecies = function(data, spec1, spec2, spec3, spec4, spec5){
  index1 = with(data, grepl(spec1, colnames(data)))
  index2 = with(data, grepl(spec2, colnames(data)))
  index3 = with(data, grepl(spec3, colnames(data)))
  index4 = with(data, grepl(spec4, colnames(data)))
  index5 = with(data, grepl(spec5, colnames(data)))
  ind = (index1 | index2 | index3 | index4 | index5) # To remove duplicated columns when there is for instance a comparaison between "m" and "A"
  species <<- makeSpeciesVector(data)
  geographies <<- makeGeographiesVector(data)
  specFac <<- selectSpeciesFac(ind) # I made sure the order of the species fac & geographies fac are correct
  geoFac <<- selectGeographiesFac(ind)
  geoNum <<- geoMap[as.character(geoFac)]
  specieColors <<- colorMap[as.character(specFac)]
  batchMatrix <<- makeBatchMatrix(data = data[ind])
  #ploidyFac <<- selectPloidiesFac(ind)
  return(data[ind])
}


## RUV (REMOVAL OF UNWANTED VARIANCE)

# ----+ function: makeRUVset
makeRUVset = function (dat){
  plotRLE(as.matrix(dat),outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  uq = betweenLaneNormalization(as.matrix(dat), which = "full")
  set = newSeqExpressionSet(uq)
  plotRLE(set,outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  
  plotPCA(set, col = specieColors, k=2, pch = geoNum, cex = 1.2, main = "First two PCs'", labels = T)
  
  return(set)
}
#---- -
## Quick analysis without any normalisation
## (Should be put in variable set)


## TODO STUFF, NORMALISATION WITH REPLICATE SAMPLE
# ----+ function: makeRUVnormalisedSet
## Normalise with genes, the "empirical" ones found with findNonDEgenes
makeRUVnormalisedSet = function(set, empircalNonDEgenes, numFactors){
  normalisedSet = RUVg(as.matrix(counts(set)), empircalNonDEgenes, k=numFactors)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, pch=geoNum, col = specieColors, main = "First two PCs'")
  return(normalisedSet)
}
#---- -
## (Should be put in variable normSet)

#Using replicate samples to remove variance
makeRUVrepNormalisedSet = function(set, empircalNonDEgenes, numFactors, replicates){
  normalisedSet = RUVs(as.matrix(counts(set)), cIdx=empircalNonDEgenes, k=numFactors, scIdx=replicates)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, pch=geoNum, col = specieColors, main = "First two PCs'")
  return(normalisedSet)
}
#---- -



# ----+ function: findNonDEgenes:
findNonDEgenes = function(model, set, numberOfGenes){
  design = model.matrix(model)
  y = DGEList(counts = counts(set), group = specFac)
  #y = calcNormFactors(y, method = "upperquartile")
  y = estimateGLMCommonDisp(y, design)
  y = estimateGLMTagwiseDisp(y, design)
  fit = glmFit(y, design)
  lrt = glmLRT(fit)
  top = topTags(lrt, n=nrow(set))$table
  #empirical <<- rownames(set)[which(!(rownames(set)%in%rownames(top)[1:numberOfGenes]))]
  return(rownames(set)[which(!(rownames(set)%in%rownames(top)[1:numberOfGenes]))]) # Should be put in variable empirical genes
  #return(rownames(top)[top$FDR >= numberOfGenes])
}

# ----+ function: findEdgeRgenes
findEdgeRgenes = function(model, group, set, pvalue){
  
  y = DGEList(counts = counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
  
  design = model.matrix(model)
  y = calcNormFactors(y, method = "upperquartile")
  y = estimateGLMCommonDisp(y, design)
  y = estimateGLMTagwiseDisp(y, design)
  plotMDS.DGEList(y, dim.plot = c(1,2), col = specieColors, cex = 1.5)
  fit = glmFit(y, design)
  lrt = glmLRT(fit, coef = 2) # Currently it's taking results of  lineal regression
  resEdgeRtopTags <<- topTags(lrt, n=nrow(set))
  resEdgeRglm <<- lrt
  de = decideTestsDGE(lrt, adjust.method="fdr", p.value = pvalue)
  de.genes = rownames(lrt)[as.logical(de)]
  set.seed(1)
  randRes = resEdgeRtopTags[rev(rownames(resEdgeRtopTags)),]
  #plotSmear(lrt, de.tags = de.genes, cex=0.5, main="MA-plot (p-value <= pvalue)")
  plot(randRes$table$logCPM, randRes$table$logFC, col=ifelse(randRes$table$FDR <= pvalue, "red", "black"), 
       cex=ifelse(randRes$table$FDR <= pvalue, 0.4, 0.3), pch = ifelse(randRes$table$FDR <= pvalue, 19, 1), xlab = "logCPM", ylab = "logFC")
  #plot(resEdgeRtopTags$table$logCPM, resEdgeRtopTags$table$logFC, col=ifelse(resEdgeRtopTags$table$FDR <= pvalue, red, black), cex=ifelse(resEdgeRtopTags$table$FDR <= pvalue, 0.4, 0.3), pch = ifelse(resEdgeRtopTags$table$FDR <= pvalue, 19, 1))
  #plotBCV(y)
  abline(h=c(-2,2), col="blue")
  abline(h=0)
  #VolcanoPlot(resEdgeRtopTags$table, 2, pvalue)
  return(resEdgeR) # return a TopTags-class object
  return(fit) # return the DGEGLM-class object
}

EdgeRresiduals = function(fit, residualType){ #residualType is "pearson","deviance" or "midp-quantile" depending on wish
  #fit is a DGEGLM-class object retured by my function findEdgeRgenes
  y = fit$counts 
  mu = fit$fitted.values
  phi = fit$dispersion 
  v = mu*(1+phi*mu)
  if(residualType == "pearson"){return((y-mu) / sqrt(v))}
  if(residualType == "deviance"){d = nbinomUnitDeviance(y,mu,phi); return(sign(y-mu) * sqrt(d))}
  if(residualType == "midp-quantile"){return(zscoreNBinom(y,mu,size=1/phi))}
}


# MAIN
resave <- function(..., list = character(), file) {
  previous  <- load(file)
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  save(list = unique(c(previous, var.names)), file = file)
}