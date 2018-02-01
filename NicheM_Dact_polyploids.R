#R script Dactylorhiza Niche model divergence
# This script downloads species coordinates from gbif. 
# it clips the pH information from SoilGrids1km
# it clip the wordclim data

library(dismo)
library('rgbif')
library(rgdal)
library(raster)
library(gtools)
library(rjson)
library(sp)
library(GSIF)
library(maps)
library(rgeos)
library(ade4)
library(caret)

###############################33
# VCF

coordinatesdact <- read.table("/home/fbalao/Datos/R/Rpackages/EnvGenExp/polyploidCoordinates.txt", header = T)
coordinates(coordinatesdact)<- ~longitude+ latitude
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(coordinatesdact) <- crs.geo 

vcf<-raster("/home/fbalao/Datos/GPS/bioclim/MODIS/out.tif")
coordinatesdact$vcf<-extract(vcf,coordinatesdact)
species<-as.factor(substr(coordinatesdact$Population, nchar(as.character(coordinatesdact$Population))-1,nchar(as.character(coordinatesdact$Population))))
boxplot(coordinatesdact$vcf~species, col=c(2,3))


#majalis data
Dfucloc<-occ_search(scientificName = c("Dactylorhiza majalis subsp. majalis","Dactylorhiza majalis"),hasGeospatialIssue=FALSE,fields='minimal', limit=50000 ,hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
gbifmap(Dfucloc[[3]])
dups <- duplicated(Dfucloc[[3]][,c(3,4)])
fuc_cleaned<- Dfucloc[[3]][!dups, ]
fuc_cleaned
coofuc<-fuc_cleaned[,c(4,3)]
# #Adding Ovidiu's coordinates
# newfuccoo<-read.table("/home/fbalao/Datos/ARTICULO/DactylorhizaDiploids/Geo/DfuchsiicoordinatesOvidiu.txt", header=F)
# colnames(newfuccoo)<- colnames(coofuc)
# rownames(newfuccoo)<-paste("i",1:14, sep="")
# #Adding Mikael's coordinates
# new2fuccoo<-read.table("/home/fbalao/Datos/ARTICULO/DactylorhizaDiploids/Geo/Mikaelloc.csv", header=T)
# 
# mikael<- within(new2fuccoo, {
#   dms <- do.call(rbind, strsplit(as.character(Lat), ":"))
#   Latitude <- as.numeric(dms[,1]) +  (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60
#   dms <- do.call(rbind, strsplit(as.character(Long.), ":"))
#   Longitude <- as.numeric(dms[,1]) +  (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60
#   })
# mikael<-mikael[,-c(3,4,7:9)]
# rownames(mikael)<-paste("mik", mikael[,1], sep="")
# mikaelfuc<-mikael[mikael$Maintaxon=="fuc",c(3,4)]
# colnames(mikaelfuc)<- colnames(coofuc)
# 
# coofuc2<-rbind(coofuc,newfuccoo, mikaelfuc)
coofuc2<-coofuc
dupsfuc <- duplicated(coofuc2)
coofuc_final<- coofuc2[!dupsfuc, ]
coordinates(coofuc_final)<- ~decimalLongitude+ decimalLatitude
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(coofuc_final) <- crs.geo 
#Remove ocurrences within 10km2
r <- raster(coofuc_final)
res(r) <- 0.008333333
r <- extend(r, extent(r)+1)
coofuc_final_sel <- as.data.frame(gridSample(coofuc_final, r, n=1))
coordinates(coofuc_final_sel)<- ~decimalLongitude+ decimalLatitude
proj4string(coofuc_final_sel) <- crs.geo 
# p <- rasterToPolygons(r)
# plot(p, border='gray')
# points(coofuc_final)
# points(acsel, cex=1, col='red', pch='x')



# traunsteineri data
Dincloc<-occ_search(scientificName = "Dactylorhiza traunsteineri",hasGeospatialIssue=FALSE,fields='minimal', limit=50000 ,hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
gbifmap(Dincloc[[3]])
dups <- duplicated(Dincloc[[3]][,c(3,4)])
inc_cleaned<- Dincloc[[3]][!dups, ]
inc_cleaned
cooinc<-inc_cleaned[,c(4,3)]

# #Adding Ovidiu's coordinates
# newinccoo<-read.table("/home/fbalao/Datos/ARTICULO/DactylorhizaDiploids/Geo/DincarnatacoordinatesOvidiu.txt", header=F)
# colnames(newinccoo)<- colnames(cooinc)
# rownames(newinccoo)<-paste("i",1:7, sep="")
# #Adding Mikael's coordinates
# 
# mikaelinc<-mikael[mikael$Maintaxon=="inc",c(3,4)]
# colnames(mikaelinc)<- colnames(cooinc)
# 
# cooinc2<-rbind(cooinc,newinccoo, mikaelinc)
cooinc2<-cooinc
dupsinc <- duplicated(cooinc2)
cooinc_final<- cooinc2[!dupsinc, ]
coordinates(cooinc_final)<- ~decimalLongitude+ decimalLatitude
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(cooinc_final) <- crs.geo 

#Remove ocurrences within 1km2
r <- raster(cooinc_final)
res(r) <- 0.008333333  #10km=0.08333333
r <- extend(r, extent(r)+1)
cooinc_final_sel <- as.data.frame(gridSample(cooinc_final, r, n=1))
coordinates(cooinc_final_sel)<- ~decimalLongitude+ decimalLatitude
proj4string(cooinc_final_sel) <- crs.geo 
# p <- rasterToPolygons(r)
# plot(p, border='gray')
# points(coofuc_final)
# points(cooinc_final_sel, cex=1, col='red', pch='x')




#MAP
layout(matrix(c(1,2), nrow=1))
map('world', xlim=c(-11,50), ylim=c(35,72))
box()
title(main="D. majalis ")
points(coofuc_final_sel, pch=16, cex=0.6, col="blue")
map('world', xlim=c(-11,50), ylim=c(35,72))
box()
points(cooinc_final_sel, pch=16, col="red", cex=0.6)
title(main="D. traunsteineri ")

#points(merged[,c(2,3)],col=c("blue","red")[merged[,1]], pch=16, cex=0.5)



#BIOCLIM DATA
#Read the raster, stack the raster, crop and merge
list.ras <- mixedsort(list.files("/home/fbalao/Datos/Dactylorhiza/BDGeoDact/wc0.5/", full.names = T, pattern = ".bil"))
bio <- stack(list.ras)
projection(bio)<-crs.geo
map('world')
newext <- drawExtent()
bio.c<-crop(bio, newext)
saveRDS(bio.c, file="/home/fbalao/Datos/Dactylorhiza/BDGeoDact/wc0.5/Bioclimvariables.rds")
#Next time I just need to load the Bioclimvariables.rds
readRDS("/home/fbalao/Datos/Dactylorhiza/BDGeoDact/wc0.5/Bioclimvariables.rds")
bio.brick <- brick(bio.c)

#Extract the bioclim variables
prevfuc<-extract(bio.c,coofuc_final_sel)
previnc<-extract(bio.c,cooinc_final_sel)
  
#Tree coverage
vcf<-raster("/home/fbalao/MODIS/out.tif")
vcffuc<-extract(vcf,coofuc_final_sel)
vcfinc<-extract(vcf,cooinc_final_sel)
vcffuc[vcffuc==200]<-0
vcfinc[vcfinc==200]<-0
boxplot(vcffuc[vcffuc<100],vcfinc[vcfinc<100], names=c("D. fuchsii", "D. incarnata"), col=c(2,3), ylab="% Tree coverage")
t.test(vcffuc[vcffuc<100],vcfinc[vcfinc<100])


biofuc<-cbind(as.data.frame(coofuc_final_sel),prevfuc,vcffuc)
bioinc<-cbind(as.data.frame(cooinc_final_sel),previnc, vcfinc)


colnames(biofuc)[22]<-'vcf'
colnames(bioinc)[22]<-'vcf'

#PCA

#1. Merge datasets
merged<-rbind(biofuc,bioinc)
sp<-c(rep('fuc',dim(biofuc)[1]), rep('inc',dim(bioinc)[1]))
merged<-cbind(sp,merged)
mergedtable<-merged
head(mergedtable)

library(zoo)
mergedtable<-na.aggregate(mergedtable)

mergedtable[] <- lapply(mergedtable, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
mergedtable
#2.Remove highly correlated variables
library('caret') 
rem<-findCorrelation(cor(mergedtable[,-c(1,2,3)]), cutoff = .80, verbose = F, names=F)
rem<-rem+3
mergedtable.c<-mergedtable[, -c(rem)]

# Background 500Km buffer polygons
coordinates(mergedtable)<- ~decimalLongitude+ decimalLatitude
proj4string(mergedtable) <- crs.geo 
mask <- raster(mergedtable)
res(mask) <- 0.008333333
x <- circles(mergedtable, d=150000, lonlat=TRUE)
pol <- gUnaryUnion(x@polygons)
samp1 <- spsample(pol, 1000, type='random', iter=25)
# get unique cells
cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy <- xyFromCell(mask, cells)
map('world', xlim=c(-11,50), ylim=c(35,72))
lines(pol, axes=TRUE, add=T, col="orchid")
points(xy, cex=0.75, pch=20, col='blue')
points(mergedtable, col='red',pch=20, cex=0.6)

backgroundclim<-extract(bio.c,xy)
backgroundvcf<-extract(vcf,xy)
backgrounddat<-cbind("background",as.data.frame(xy),backgroundclim, backgroundvcf)
backgrounddat<-na.omit(backgrounddat)
backgrounddat.c<-backgrounddat[, -c(rem)]
colnames(backgrounddat.c)<-colnames(mergedtable.c)

# Background for each species 100Km buffer polygons
#incarnata background
mask <- raster(cooinc_final_sel)
res(mask) <- 0.008333333
x <- circles(cooinc_final_sel, d=100000, lonlat=TRUE)
polinc <- gUnaryUnion(x@polygons)
samp1 <- spsample(polinc, 3000, type='random', iter=25)
# get unique cells
cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)
incxy <- xyFromCell(mask, cells)
plot(polinc, axes=TRUE, add=T)
points(xy, cex=0.75, pch=20, col='blue')
points(mergedtable[sp=="inc",], col='red',pch=20, cex=0.6)


incbackgroundclim<-extract(bio.c,incxy)
incbackgroundvcf<-extract(vcf,incxy)
incbackgrounddat<-cbind("incbackground",as.data.frame(incxy),incbackgroundclim, incbackgroundvcf)
incbackgrounddat<-na.omit(incbackgrounddat)
incbackgrounddat.c<-incbackgrounddat[, -c(rem)]
colnames(incbackgrounddat.c)<-colnames(mergedtable.c)



#fuchsii background
mask <- raster(coofuc_final_sel)
res(mask) <- 0.008333333
x <- circles(coofuc_final_sel, d=100000, lonlat=TRUE)
polfuc <- gUnaryUnion(x@polygons)
samp1 <- spsample(polfuc, 3000, type='random', iter=25)
# get unique cells
cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)
fucxy <- xyFromCell(mask, cells)
lines(polfuc, axes=TRUE, col=2)
points(fucxy, cex=0.75, pch=20, col='blue')
points(coofuc_final_sel, col='red',pch=20, cex=0.6)
points(cooinc_final_sel, col='green',pch=20, cex=0.6)

fucbackgroundclim<-extract(bio.c,fucxy)
fucbackgroundvcf<-extract(vcf,fucxy)
fucbackgrounddat<-cbind("fucbackground",as.data.frame(fucxy),fucbackgroundclim, fucbackgroundvcf)
fucbackgrounddat<-na.omit(fucbackgrounddat)
fucbackgrounddat.c<-fucbackgrounddat[, -c(rem)]
colnames(fucbackgrounddat.c)<-colnames(mergedtable.c)


plot(pol, axes=TRUE)
points(fucxy, cex=0.75, pch=20, col='blue')
points(incxy, cex=0.75, pch=20, col='grey55')
points(coofuc_final_sel, col='red',pch=20, cex=0.6)
points(cooinc_final_sel, col='green',pch=20, cex=0.6)





#Merge background and species data

todo<-rbind(mergedtable.c,fucbackgrounddat.c, incbackgrounddat.c)

#Niche overlapping

#number of interation for the tests of equivalency and similarity
iterations<-100
#resolution of the gridding of the climate space
R=100
#################################### PCA-ENVIRONMENT ##################################
data<-todo[,-c(1:3)]
library(ade4)
w<-c(rep(0,nrow(mergedtable.c)),rep(1,nrow(fucbackgrounddat.c)), rep(1,nrow(incbackgrounddat.c)))
pca.cal <-dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
barplot(pca.cal$eig)
gcol = c("blue", "red", rgb(0,1,0,0.3), rgb(1,0.5,0,0.3))

scatter(pca.cal, clab.row = 0, posieig = "none", col=NULL)
s.class(pca.cal$li, todo[,1], col = gcol, add.plot = TRUE, cstar = 0, clabel = 0, cellipse = 1.5)


  row.fuc<-which(todo[,1] == "fuc") # rows in data corresponding to sp1
  row.inc<-which(todo[,1] == "inc")# rows in data corresponding to sp2
  row.bacfuc<-which(todo[,1] == "fucbackground") # rows in data corresponding to sp1
  row.bacinc<-which(todo[,1] == "incbackground")# rows in data corresponding to sp2
                                               
                 
  # predict the scores on the axes
  scores.clim<- pca.cal$li[(nrow(mergedtable.c)+1):nrow(data),]  #scores for global climate
  scores.fuc<- pca.cal$li[row.fuc,]					#scores for sp1
  scores.inc<- pca.cal$li[row.inc,]					#scores for sp2
  scores.bacfuc<- pca.cal$li[row.bacfuc,]					#scores for sp1
  scores.bacinc<- pca.cal$li[row.bacinc,]	
  

library(ecospat)
zfuc<- ecospat.grid.clim.dyn(scores.clim, scores.bacfuc, scores.fuc,R=100)
zinc<- ecospat.grid.clim.dyn(scores.clim, scores.bacinc, scores.inc,R=100)
#POR AQUI!!!!
similarity.testfucinc<-ecospat.niche.similarity.test (zfuc, zinc, 100, one.sided=F)
similarity.testincfuc<-ecospat.niche.similarity.test (zinc, zfuc, 100, one.sided=F)
equivalency.test<-ecospat.niche.equivalency.test (zfuc, zinc, 100)

ecospat.plot.overlap.test(similarity.testincfuc, type="D", title="Similarity test D. inc -> D. fuc")
ecospat.plot.overlap.test(similarity.testfucinc, type="D", title="Similarity test D. fuc -> D. inc")

ecospat.plot.niche.dyn(  zfuc, zinc, title="Niche Categories and Species Density",quant=0.75)

ecospat.fun.arrows (scores.fuc, scores.inc, scores.clim, scores.clim)

#Test divergence
PC1dn<-t.test(scores.fuc[,1],scores.inc[,1])
PC2dn<-t.test(scores.fuc[,2], scores.inc[,2])

dbPC1<-vector()
for (i in 1:9999){
  dbiPC1<-mean(sample(scores.bacinc[,1],500, replace=T))
  dbfPC1<-mean(sample(scores.bacfuc[,1],500, replace=T))
  dbPC1[i]<-dbfPC1-dbiPC1
  }

dbPC2<-vector()
for (i in 1:9999){
  dbiPC2<-mean(sample(scores.bacinc[,2],500, replace=T))
  dbfPC2<-mean(sample(scores.bacfuc[,2],500, replace=T))
dbPC2[i]<-dbfPC2-dbiPC2
  }

dn1<-mean(scores.fuc[,1])- mean(scores.inc[,1])
db1<- quantile(dbPC1,c(0.025, .975))

dn2<-mean(scores.fuc[,2])- mean(scores.inc[,2])
db2<-quantile(dbPC2,c(0.025, .975))

#dbIC= (-0.011362342  0.002173957)
#daPC1=0.256439
#3.PCA on the climatic variables
library(ade4)
pca1 <- dudi.pca(mergedtable.c[,-c(1,2,3)], scann = FALSE, nf = 3)
barplot(pca1$eig)
gcol = c("blue", "red")

s.class(dfxy = pca1$li, fac = mergedtable[,1], col = gcol, xax = 1, yax = 2)
s.corcircle(pca1$co, xax = 1, yax = 2)

scatter(pca1, clab.row = 0, posieig = "none")
s.class(pca1$li, mergedtable[,1], col = gcol, add.plot = TRUE, cstar = 0, clabel = 0, cellipse = 1.5)

#MANOVA
manovabio<-manova(as.matrix(mergedtable.c[,-c(1,2,3)]) ~ as.factor(mergedtable.c[,1]) )
summary(manovabio)

#Niche overlaping
library(ecospat)

#Independent ANOVA
a1<-aov(mergedtable.c$bio_1 ~ as.factor(mergedtable.c[,1]))
a2<-aov(mergedtable.c$bio_2 ~ as.factor(mergedtable.c[,1]))
a3<-aov(mergedtable.c$bio_3 ~ as.factor(mergedtable.c[,1]))
a4<-aov(mergedtable.c$bio_8 ~ as.factor(mergedtable.c[,1]))
a5<-aov(mergedtable.c$bio_10 ~ as.factor(mergedtable.c[,1]))
a6<-aov(mergedtable.c$bio_15 ~ as.factor(mergedtable.c[,1]))
a7<-aov(mergedtable.c$bio_18 ~ as.factor(mergedtable.c[,1]))

anova(a1,a2,a3)
library("xtable")
print(xtable(anova(a1), type="html",file="anova_bio1.html"))
print(xtable(anova(a2), type="html",file="anova_bio2.html"))
print(xtable(anova(a3), type="html",file="anova_bio3.html"))
print(xtable(anova(a4), type="html",file="anova_bio8.html"))
print(xtable(anova(a5), type="html",file="anova_bio10.html"))
print(xtable(anova(a6), type="html",file="anova_bio15.html"))
print(xtable(anova(a7), type="html",file="anova_bio18.html"))

layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2, byrow =T ))
boxplot(mergedtable.c$bio_1 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio1")
boxplot(mergedtable.c$bio_2 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio2")
boxplot(mergedtable.c$bio_3 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio3")
boxplot(mergedtable.c$bio_8 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio8")
boxplot(mergedtable.c$bio_10 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio10")
boxplot(mergedtable.c$bio_15 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio15")
boxplot(mergedtable.c$bio_18 ~ as.factor(mergedtable.c[,1]),col = gcol, main="bio18")
plot.new()

#MANOVA in PC
manovabiopca<-manova(cbind(pca1$li$Axis1,pca1$li$Axis2, pca1$li$Axis3) ~ as.factor(mergedtable.c[,1]) )
summary(manovabiopca)

layout(matrix(c(1,1,1,2,3,4),nrow=2, byrow =T ))
scatter(pca1, clab.row = 0, posieig = "none")
s.class(pca1$li, mergedtable[,1], col = gcol, add.plot = TRUE, cstar = 0, clabel = 0, cellipse = 1.5)
boxplot(pca1$li$Axis1~as.factor(mergedtable.c[,1]), col=gcol, ylab="Component 1",xlab="Species")
boxplot(pca1$li$Axis2~as.factor(mergedtable.c[,1]), col=gcol, ylab="Component 2",xlab="Species")
boxplot(pca1$li$Axis3~as.factor(mergedtable.c[,1]), col=gcol,  ylab="Component 3",xlab="Species")


#Soil
soilgrids.r <- REST.SoilGrids(c("ORCDRC","PHIHOX"))
ovfuc <- over(soilgrids.r, coofuc_final_sel)
ovinc <- over(soilgrids.r, cooinc_final_sel)
#pH
library(lattice)
library(aqp)
data(soil.legends)
plot(density(na.omit(ovfuc$PHIHOX.M.sd1)/10), col=3)
lines(density(na.omit(ovinc$PHIHOX.M.sd1)/10), col=2)

PHIHOX.range = range(soil.legends[["PHIHOX"]]$MIN, soil.legends[["PHIHOX"]]$MAX)

