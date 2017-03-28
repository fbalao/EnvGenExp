compareEdgeRresults<-function(var){
res0.1<-read.table(paste0("./rescpm0.1/",var,".txt"), header=T,sep=",", row.names = 1)
res1<-read.table(paste0("./respcpm1/",var,".txt"), header=T,sep=",", row.names = 1)
res2<-read.table(paste0("./respcpm2//",var,".txt"), header=T,sep=",", row.names = 1)
sigres0.1<-rownames(res0.1[res0.1[,5]<0.05,])
sigres1<-rownames(res1[res1[,5]<0.05,])
sigres2<-rownames(res2[res2[,5]<0.05,])
nsigres0.1<-length(sigres0.1)
nsigres1<-length(sigres1)
nsigres2<-length(sigres2)

c0.1_1<-Reduce(intersect, list(sigres0.1,sigres1))
c0.1_2<-Reduce(intersect, list(sigres0.1,sigres2))
c1_2<-Reduce(intersect, list(sigres1,sigres2))
common<-Reduce(intersect, list(sigres0.1,sigres1,sigres2))
ncommon<-length(common)
nc0.1_1<-length(c0.1_1)
nc0.1_2<-length(c0.1_2)
nc1_2<-length(c0.1_2)

require(VennDiagram)

venn.plot <- draw.triple.venn(
  area1 = nsigres0.1,
  area2 = nsigres1,
  area3 = nsigres2,
  n12 = nc0.1_1,
  n23 = nc1_2,
  n13 = nc0.1_2,
  n123 = ncommon,
  category = c("CPM<0.1", "CPM<1", "CPM<2"),
  fill = c("blue", "red", "green"),
  scaled=TRUE)

tiff(filename = paste0(var,"_Venn_diagram.tiff"), compression = "lzw")
grid.draw(venn.plot);
dev.off();
}

compareEdgeRresults("bio_3")
compareEdgeRresults("bio_13")
compareEdgeRresults("bio_15")
compareEdgeRresults("bio_5")
compareEdgeRresults("pH")
barplot(c(nsigres0.1,nsigres1,nsigres2), names.arg = c("cpm<0.1","cpm<1","cpm<2"))
