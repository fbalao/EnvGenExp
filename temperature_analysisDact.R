library(data.table)
library(ggplot2)
library(xts)
library(edgeR)

setwd("/home/fbalao/Datos/Dactylorhiza/environment/temperature")

#Take names of dataloggers and populations
filenames <- list.files(path=".",
                        pattern="*.csv", recursive = T)
tmp<-strsplit(filenames, split='/')
loggerinfo<-do.call(rbind, lapply(tmp, rbind)) #a dataframe with dataloggers names and populations
loggerinfo[,3]<-substr(loggerinfo[,3],1,nchar(loggerinfo[,3])-4)

#Read all the csv
sources.files  <- list.files(path='.',
                               recursive=T,
                               pattern='*.csv'
                               ,full.names=T)
temp <- lapply(sources.files, fread, sep=",", skip=1L)
names(temp)<-paste(loggerinfo[,2], loggerinfo[,3], sep="_")

datos <- rbindlist(temp, idcol=TRUE)
datos<-datos[,c(1,3,4)]
colnames(datos)<-c("ID","Date","Temperature")


datos$Date<-gsub('\\.', '-', datos$Date)
datos$Date<-as.POSIXct(datos$Date,tz = "GMT", format ="%d-%m-%Y %H:%M:%OS")

# colnames(datos)[2]<-"Date"
# plot(datos$Date,datos$Temp,xaxt="n",ylab="ºC", type="l") 
# daterange<-c(min(datos$Date),max(datos$Date))
# axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="month"), format="%b")

datos$species <- substr(datos$ID,1,3) 
datos$pop <- substr(datos$ID,5,9) 


#Reorder the data.frame
datos2<-datos[,c(2,3,4,5,1)]

# Bioclim estimates
bio1<- function(x){mean(apply.monthly(x, FUN=mean))}
bio2<- function(x){
monthMax<-apply.monthly(x, FUN=max)
monthMin<-apply.monthly(x, FUN=min)
mean(monthMax-monthMin)
}
bio4<-function(x){sd(apply.monthly(x, FUN=mean))}
bio5<-function(x){max(apply.monthly(x, FUN=max))}
bio6<-function(x){min(apply.monthly(x, FUN=min))}
bio7<-function(x){bio5(x)-bio6(x)}
bio3 <-function(x){(bio2(x)*100)/bio7(x)}
bio10<-function(x){min(apply.quarterly(x, FUN=mean))}
bio11<-function(x){max(apply.quarterly(x, FUN=mean))}


#Loop to create xts objects for each datalogger

data<-split(datos2, by="ID")
res<-data.frame()
for (i in 1: length(data)){
subdata<-as.xts(data[[i]][,c(1,2)])
subdata["/2015 & /2016"] # Comment if you want the complet dataset
res[i,1]<-bio1(subdata)
res[i,2]<-bio2(subdata)
res[i,3]<-bio3(subdata)
res[i,4]<-bio4(subdata)
res[i,5]<-bio5(subdata)
res[i,6]<-bio6(subdata)
res[i,7]<-bio7(subdata)
res[i,8]<-bio10(subdata)
res[i,9]<-bio11(subdata)
}

colnames(res)<-c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio10","bio11")
res$ID<-names(data)
res$species <- substr(res$ID,1,3) 
res$pop <- substr(res$ID,5,9) 

meanres<-with(res, aggregate(cbind(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio10,bio11)~as.factor(res$species):as.factor(res$pop), FUN=mean))


####################################3
#Loop to create xts objects for each datalogger

# data<-split(datos2, by="ID")
# res<-data.frame()
# for (i in 1: length(data)){
#   subdata<-as.xts(data[[i]][,c(1,2)])
#   subdata<-subdata[.indexmon(subdata)%in% c(5)]
#  res[i,1]<-mean(apply.monthly(subdata,FUN = max))
#  res[i,2]<-mean(apply.monthly(subdata,FUN = min))}
# 
# colnames(res)<-c("maxJune","minJune")
# res$ID<-names(data)
# res$species <- substr(res$ID,1,3) 
# res$pop <- substr(res$ID,5,9) 
# 
# meanres<-with(res, aggregate(cbind(maxJune,minJune)~as.factor(res$species):as.factor(res$pop), FUN=mean))
# 

#######################333333
# Some analysis and plots


diftes<-numeric()
for (i in 1:9){
tetst<-t.test(res[,i]~res[,11])
diftes[i]<-tetst$p.value
}

boxplot(res$bio4~res$pop, col=c(2,3)[as.factor(res$species)])



###################################
### SOIL analysis
###############################
setwd("/home/fbalao/Datos/R/Rpackages/EnvGenExp/")
soil<-read.table("fullsoil.csv", header=T)

soil$indexsoil<-paste(substr(soil$Label_code,1,nchar(as.character(soil$Label_code))-1), soil$species, sep="_")

soilmean<-aggregate(.~indexsoil, data=soil[,-c(1,2)], FUN = function(x) mean(x, na.rm=TRUE),na.action = na.pass)
soilmean

meanres2<-meanres[-c(6,7),]

levels(meanres2$`as.factor(res$species)`)[levels(meanres2$`as.factor(res$species)`)=="maj"] <- "majalis"
levels(meanres2$`as.factor(res$species)`)[levels(meanres2$`as.factor(res$species)`)=="tra"] <- "trau"
levels(meanres2$`as.factor(res$pop)`) <- gsub("_", "", levels(meanres2$`as.factor(res$pop)`))
rownames(meanres2)<-paste(meanres2$`as.factor(res$pop)`, meanres2$`as.factor(res$species)`, sep="_")
rownames(soilmean)<-soilmean$indexsoil
## BRIT 2 stuf
soilmean2 <- rbind(soilmean, BRI_2_trau = soilmean[soilmean$indexsoil=="BRI_1_trau",])

## MERGE SOIL AND TEMP
soil.temp <- cbind(soilmean2[order(rownames(soilmean2)),], meanres2[order(rownames(meanres2)),])
rownames(soil.temp)->soil.temp$indexsoil
######################################################################
############ IMPORT COORDNATES ALLOWS MAP IND-POP LINK
#######################################################################
coordinates <- read.table("/home/fbalao/Datos/R/Rpackages/EnvGenExp/polyploidCoordinates.txt", header = T)
keycoo<-unique(coordinates$Population)

counts<-read.table("/home/fbalao/Datos/R/Rpackages/EnvGenExp/datFilteredCounts.txt")

soil.temp<-soil.temp[order(rownames(soil.temp)),]

#PH

ph<-read.table("/home/fbalao/Datos/R/Rpackages/EnvGenExp/pH.txt", header=T)
phmeans<-colMeans(ph, na.rm = T)
phmeans2<-phmeans[order(names(phmeans))]

na.omit(merge(phmeans2,soil.temp,by="row.names",all.x=TRUE))

toremove_env<-setdiff(rownames(soil.temp), levels(coordinates$Population))



envmatrix<-(merge(phmeans2,soil.temp,by="row.names",all.x=TRUE))
rownames(envmatrix)<-envmatrix$Row.names
finalenvmatrix<-envmatrix[as.character(keycoo),]
colnames(finalenvmatrix)[26]<-"species"
colnames(finalenvmatrix)[27]<-"pop"
colnames(finalenvmatrix)[2]<-"pHfield"

finalenvmatrix2<-finalenvmatrix[,-c(1,3,26,27)]
finalenvmatrix3<-finalenvmatrix2[, ! colnames(finalenvmatrix2) %in% c("Al", "Pb","Ca","Cd","Cu","Mo", "Fe", "Mg", "Mn", "Ni", "P", "S", "Zn", "Carbon."),]
# Asses multicollinearity
mtemp<-finalenvmatrix[,-c(1,3,15,19,27,26)] # Removing factor but also Mo and Pb
library('usdm')
x<-vifstep(finalenvmatrix3, th=10)

# selected<-setdiff(x@variables,x@excluded)
# reducedmatrix<-subset(mtemp, select = selected)

# Duplicate rows
ind<-table(coordinates$Population) #freq ind per pop
indsorted<-ind[rownames(finalenvmatrix3)]
finalenvmatrixrep<-finalenvmatrix3[rep(seq_len(nrow(finalenvmatrix3)), times=indsorted),]
saveRDS(finalenvmatrixrep, "finalenvmatrixrep.RDS")

##################################################################################################################3333
####################### CORRELATIONS ENV-GENE EXPRESSION 
######################
##################################################################################################################3

#load counts
filteredcounts = read.table("datFilteredCounts.txt", header = TRUE, row.names = 1) 
load("normSet.RData")
colnames(normSet$W)<-c("W1","W2","W3")
#Make the set allowing normalisation + normalisation
#set <- makeRUVset(dat = dat.filtered)
#normSet <- makeRUVrepNormalisedSet(set, rownames(counts(set)), 3, batchMatrix)

reslist<-list()

for (i in 1:dim(finalenvmatrixrep)[2]){
design <- model.matrix(~finalenvmatrixrep[,i]+normSet$W)
y <- DGEList(counts = filteredcounts, group = NULL) #Don't forget, model is taking into account the RUV normalised matrix of counts
#get an idea of how to filter
y <- calcNormFactors(y)
y <- estimateDisp(y, design = design, robust = T)
#plotBCV(y)
fit <- glmFit(y, design, robust = T)
hist(fit$coefficients[,2], breaks = 100)
res<-glmLRT(fit, coef=2, contrast=NULL)
de = decideTestsDGE(res, adjust.method="fdr", p.value = 0.05)
de.genes = rownames(res)[as.logical(de)]
#plotSmear(res, de.tags = de.genes, cex=0.5, main="MA-plot (fdr <= 0.05)")
reslist[[i]] <- topTags(res, n=nrow(filteredcounts))
}
names(reslist)<-colnames(finalenvmatrixrep)

#loadRSD("env_gen_analysis.RDS")
plot(finalenvmatrixrep$Nitrogen., normSet$normalizedCounts['comp3606_c0_seq1',] , cex=1, pch=1, ylab="Normalized counts", xlab="Nitrogen", col=1, main="comp3606_c0_seq1")


