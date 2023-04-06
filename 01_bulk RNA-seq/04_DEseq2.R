####PCA#####
library(DESeq2)
data.path="02_bulkRNA/00_count/"
sample.path="02_bulkRNA/00_all/"
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list

HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-8)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
tt=lapply(data.list, HTseq.handling,46)
head(tt[[1]])
tail(tt[[1]])
head(tt[[9]])
data.name=unique(substr(data.list,46,nchar(data.list[1])-10))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)
write.csv(data.all,paste(sample.path,'merge_sample.csv',sep = ''),row.names =F)

data.combn=data.frame(c(data.name[1],data.name[2],data.name[3]))
all.name=as.character(data.combn[,1])
all.sample=c(paste(all.name[1],c(1,2,3),sep='-'),paste(all.name[2],c(1,2,3,4,5,6,7),sep='-'),
             paste(all.name[3],c(1,2,3),sep='-'))
samplecondition=factor(c(rep(all.name[1],3),rep(all.name[2],7),rep(all.name[3],3)),
                       levels = c(all.name[1],all.name[2],all.name[3]))
samplecondition ##差异比较矩阵,说明是对照还是处理
mycol.data=data.frame(row.names = all.sample,condition=samplecondition)
mycol.data ##样品信息矩阵，即condition
mycount.data=data.all
rownames(mycount.data)=mycount.data[,1]
mycount.data=mycount.data[,-1] ##表达矩阵（所有样品）
class(mycount.data)

##构建dds矩阵
mydds=DESeqDataSetFromMatrix(countData = mycount.data,colData = mycol.data,design = ~ condition)
head(mycount.data)
head(mycol.data)
samplecondition

##标准化
mydds=DESeq(mydds)
myres=results(mydds)
########rlog
rld=rlog(mydds)
write.csv(assay(rld),paste(sample.path,"merge_sample_rlog.csv",sep=''),row.names = T)
myres=myres[order(myres$padj),]
##画图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
print('principal components analysis and plot')
pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
write.csv(pca_data,paste(sample.path,"PCA_data.csv",sep=""),row.names = F)
pca_data <- read.csv(paste(sample.path,"PCA_data.csv",sep=""))
percentVar <- round(100 * attr(pca_data, "percentVar"))
pdf(paste(sample.path,"PCA_plot.pdf",sep=""),width=4, height=3)
p=ggplot(pca_data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) + 
  scale_color_manual(values = c('#a0c3d5', '#c9caca', '#4b2415'))+ 
  scale_shape_manual(values = c(16,16,16))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+theme_classic()
print(p)
dev.off()

####组间差异基因计算#####
library(DESeq2)
all_symbol=read.csv("Homo_sapiens.GRCh37.87_chr_with_len.tsv",sep='\t')
data.path="02_bulkRNA/00_count/"
sample.path="02_bulkRNA/01_DEGs/03_F_A/"
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list
HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-8)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
tt=lapply(data.list[c(1:3,11:13)], HTseq.handling,46)
head(tt[[1]])
tail(tt[[1]])
head(tt[[6]])
data.name=unique(substr(data.list,46,nchar(data.list[1])-10))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)
data.all=data.all[,c(1,5:11,2:4)]##换列 对照在前 实验在后
write.csv(data.all,paste(sample.path,'merge_sample.csv',sep = ''),row.names =F)

data.combn=data.frame(c(data.name[1],data.name[3]))
all.name=as.character(data.combn[,1])
all.sample=c(paste(all.name[1],c(1,2,3),sep='-'),paste(all.name[2],c(1,2,3),sep='-'))
samplecondition=factor(c(rep(all.name[1],3),rep(all.name[2],3)),
                       levels = c(all.name[1],all.name[2]))
samplecondition ##差异比较矩阵,说明是对照还是处理
mycol.data=data.frame(row.names = all.sample,condition=samplecondition)
mycol.data ##样品信息矩阵，即condition
mycount.data=data.all
rownames(mycount.data)=mycount.data[,1]
mycount.data=mycount.data[,-1] ##表达矩阵（所有样品）
class(mycount.data)

##构建dds矩阵
mydds=DESeqDataSetFromMatrix(countData = mycount.data,colData = mycol.data,design = ~ condition)
head(mycount.data)
head(mycol.data)
samplecondition

##标准化
mydds=DESeq(mydds)
myres=results(mydds)
########rlog
rld=rlog(mydds)
write.csv(assay(rld),paste(sample.path,"sample_rlog.csv",sep=''),row.names = T)
myres=myres[order(myres$padj),]

##提取差异分析结果
all_gene.data=merge(as.data.frame(myres),as.data.frame(counts(mydds,normalize=TRUE)),by="row.names",sort=FALSE)
resdata=all_gene.data
colnames(all_symbol)[1]="gene_id"
colnames(resdata)[1]="gene_id"
resdata_anno=merge(all_symbol,resdata,by="gene_id")
write.csv(resdata_anno,paste(sample.path,"all_gene.csv",sep=''),row.names = F)
######diff_gene
diff_gene=subset(resdata_anno,padj<0.05&(log2FoldChange>=1.5|log2FoldChange<=-1.5))
write.csv(diff_gene,paste(sample.path,"diff_gene.csv",sep=''),row.names = F)
#####up_gene
up_gene=subset(diff_gene,log2FoldChange>=0.25)
write.csv(up_gene,paste(sample.path,"up_gene.csv",sep=''),row.names = F)
#####down_gene
down_gene=subset(diff_gene,log2FoldChange<=-0.25)
write.csv(down_gene,paste(sample.path,'down_gene.csv',sep=''),row.names = F)
