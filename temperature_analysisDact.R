library(data.table)
library(ggplot2)
library(xts)

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