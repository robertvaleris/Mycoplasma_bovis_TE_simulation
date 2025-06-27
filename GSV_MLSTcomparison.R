data<-read.csv(file="/Users/robertjvalerischacin/Downloads/supplementary_materials_new.csv", stringsAsFactors = F, header = T)

data2<-data[1:620,]

data2<-data2[!grepl("closest", data2$MLST),]

data2<-data2[!grepl("Closest", data2$MLST),]

data2<-data2[!is.na(data2$MLST),]

data2$MLST<-as.numeric(data2$MLST)

library(igraph)
library(mcclust)


igraph::compare(data2$GSV, data2$MLST, method="nmi")
igraph::compare(data2$GSV, data2$MLST, method="rand")
igraph::compare(data2$GSV, data2$MLST, method="adjusted.rand")
vi.dist(data2$GSV, data2$MLST)