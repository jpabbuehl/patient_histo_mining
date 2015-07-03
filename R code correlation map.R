```{r}
HighTumorInfLN_data<-read.table('High tumor infiltrated LN.csv',head=TRUE,sep=",",row.names=NULL,stringsAsFactors = FALSE)
id<-c(9:ncol(HighTumorInfLN_data))
HighTumorInfLN_data[,id]<-as.numeric(as.character(unlist(HighTumorInfLN_data[,id])))

LowTumorInfLN_data<-read.table('Low tumor infiltrated LN.csv',head=TRUE,sep=",",row.names=NULL,stringsAsFactors = FALSE)
id<-c(8:ncol(LowTumorInfLN_data))
LowTumorInfLN_data[,id]<-as.numeric(as.character(unlist(LowTumorInfLN_data[,id])))

LN_responder<-read.table('LN from Responder.csv',head=TRUE,sep=",",row.names=NULL,stringsAsFactors = FALSE)
id<-c(8:ncol(LN_responder))
LN_responder[,id]<-as.numeric(as.character(unlist(LN_responder[,id])))

LN_nonresponder<-read.table('LN from Non Responder.csv',head=TRUE,sep=",",row.names=NULL,stringsAsFactors = FALSE)
id<-c(9:ncol(LN_nonresponder))
LN_nonresponder[,id]<-as.numeric(as.character(unlist(LN_nonresponder[,id])))

data <- lapply(full_data[,9:ncol(full_data)], FUN = function(x) as.numeric(gsub("%", "", x)))
data<- as.data.frame(lapply(full_data[,8:ncol(full_data)],function(x) if(is.character(x)|is.factor(x)) gsub("%","",x) else x))
data<-as.matrix(data)
class(data)<-'numeric'

nb_na<-apply(data,2,function(x) sum(is.na(x)))
subset<-data[,nb_na<30]

library(Hmisc)
library(corrplot)
correlation<-rcorr(subset)

cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

confidence<-cor.mtest(correlation$r, 0.95)[[1]]
rownames(confidence)<-rownames(correlation)
colnames(confidence)<-colnames(correlation)

col1 <- colorRampPalette(c("green", "white", "red"))
tiff(filename = "final.tiff", width =2000, height = 2000,
     units = "px", pointsize = 12, bg = "white", res = NA)
corrplot(correlation$r,method="color",
         type="full",
         p.mat = confidence, insig = "blank",order="hclust", col=col1(100),addrect=11,
         tl.col = "black",addgrid.col = "gray50",diag=FALSE,tl.cex=1, cl.pos = "b",cl.cex=3,cl.offset=0.1)
dev.off()

col1 <- colorRampPalette(c("green", "white", "red"))
tiff(filename = "final without pval.tiff", width =2000, height = 2000,
     units = "px", pointsize = 12, bg = "white", res = NA)
corrplot(correlation$r,method="color",
         type="full",insig = "blank",order="hclust", col=col1(100),addrect=11,
         tl.col = "black",addgrid.col = "gray50",diag=FALSE,tl.cex=1, cl.pos = "b",cl.cex=3,cl.offset=0.1)
dev.off()

col1 <- colorRampPalette(c("green", "black", "red"))
tiff(filename = "final second color.tiff", width =2000, height = 2000,
     units = "px", pointsize = 12, bg = "white", res = NA)
corrplot(correlation$r,method="color",
         type="full",
         p.mat = confidence, insig = "blank",order="hclust", col=col1(100),addrect=11,
         tl.col = "black",addgrid.col = "gray50",diag=FALSE,tl.cex=1, cl.pos = "b",cl.cex=3,cl.offset=0.1)
dev.off()

col1 <- colorRampPalette(c("green", "black", "red"))
tiff(filename = "final second color without pval.tiff", width =2000, height = 2000,
     units = "px", pointsize = 12, bg = "white", res = NA)
corrplot(correlation$r,method="color",
         type="full",insig = "blank",order="hclust", col=col1(100),addrect=11,
         tl.col = "black",addgrid.col = "gray50",diag=FALSE,tl.cex=1, cl.pos = "b",cl.cex=3,cl.offset=0.1)
dev.off()