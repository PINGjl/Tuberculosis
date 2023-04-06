library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)

seu.sub <- readRDS(file = "/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
deg = read.csv('/05_Result/07_co_lung/01_ref/06_IF/dis_time_TF/up_down_module_degs.csv')

for (i in c("AT1","AT2","Ciliated","Clu_Gob", "EC","Stroma",'Mast',"Classical_MC","Nonclassical_MC","Mac1","Mac2","Neu","mDC","pDC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')){
tmp <- subset(seu.sub, celltype==i)
exp.mat <- GetAssayData(tmp, slot='data')
deg_tmp=subset(deg, celltype==i)
genes <- unique(as.character(deg_tmp$gene))
exp.mat <- as.matrix(exp.mat)
exp.mat <- as.matrix(exp.mat[genes,])

setwd(paste0("/05_Result/07_co_lung/01_ref/06_IF/dis_time_TF/",i,"/"))

org="hgnc"
dbDir="/03_Database/02_SCENIC/hg19_cisTarget_databases"
myDatasetTitle="SCENIC analysis of dis_time"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores线程

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# Gene filter
genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exp.mat),
                           minSamples=ncol(exp.mat)*.01)
exprMat_filtered <- exp.mat[genesKept, ]

# Run Genie3
runCorrelation(exprMat_filtered, scenicOptions)
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
        print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------

# Run the remaining
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 15
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions)
        runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        print(paste0('----------------------',i,'_data completed','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------
}
