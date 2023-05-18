# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Transcriptomic features associated with response
# Section:      Results - Tumour proliferation and immune signatures
#=======================================================================================

rm (list=ls())

#load packages
library (data.table)
library (edgeR)
library (EnsDb.Hsapiens.v86)
library (fgsea)
library (ggbiplot)
library (ggplot2)
library (ggpmisc)
library (ggridges)
library (MASS)
library (org.Hs.eg.db)
library (Pigengene)
library (ReactomePA)
library (readxl)
library (sm)
library (stringr)
library (viridis)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "D:/final project/neoadjuvant-therapy-response-predictor-master/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

figure_font_size=12

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status <- ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status <- ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)
metadataFull <- metadata

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]

# load list of breast cancer driver genes from resources directory
driverGenes<- scan(paste0(resourcesDir,"breast-cancer-driver-genes.txt"), what=character(),skip = 1)

# load Gene Ensembl ID to Hugo ID dictionary
ensemblToHugo <- read.table(paste0(resourcesDir,"EnsemblID.to.Hugo.v87.tsv.gz"), header=T, stringsAsFactors = F,sep="\t")

# load RNA data (Supplementary Table 3)
rnadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 3))


#=========================================================================
# Differential gene expression and enrichment
# Figure 3a
#=========================================================================

p <- metadata

# load RNAseq raw counts (Methods)
transneo.counts <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-RNAseq-rawcounts.tsv.gz"),header=T, sep="\t",stringsAsFactors = F),row.names = 1)
transneo.counts <- transneo.counts[,colnames(transneo.counts) %in% p$Donor.ID]

# update metadata - retain only samples that have RNAseq data
p <- p[p$Donor.ID %in% colnames(transneo.counts),]

# sanity check
stopifnot(sum(p$Donor.ID!=colnames(transneo.counts))==0)

# 149 tumours have RNAseq data, and associated RCB assessment + had adequate chemotherapy exposure
dim(p)

# construct differential expression matrix
y <- DGEList(transneo.counts)
minCPM      <- 1
minNoToKeep <- 10
keep        <- rowSums(cpm(y)>minCPM)>=minNoToKeep
y <- y[keep , , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")

# construct design matrix
rcb    <- factor(p$pCR.RD, levels=c("RD","pCR"))
design <- model.matrix(~rcb, data = y$samples)
coef   <- ncol(design)
head(design,2)

# perform DE: pCR vs RD - as per edgeR documentation
y       <- estimateDisp(y, design = design, robust = TRUE)
fit     <- glmQLFit(y, design, robust=TRUE)
qlf     <- glmQLFTest(fit, coef=coef)
results <- as.data.frame(topTags(qlf,n = Inf))
is.de   <- decideTestsDGE(qlf, p.value=0.05)

# 2,071 genes are under-expressed, and 2,439 genes are over-expressed in tumours attaining pCR (FDR<0.05). 
summary(is.de)

# Annotate full gene list with their Hugo Gene name
annotatedResults <- results
annotatedResults[rownames(annotatedResults) %in% rownames(y)[is.de[,1]==1],"expression"] <- "overexpressed"
annotatedResults[rownames(annotatedResults) %in% rownames(y$counts)[is.de[,1]==-1],"expression"] <- "underexpressed"
annotatedResults[rownames(annotatedResults) %in% rownames(y$counts)[is.de[,1]==0],"expression"] <- "notDE"
annotatedResults <- merge(x=annotatedResults,y=ensemblToHugo, by.x=0,by.y=1, sort=F)


cut_off_fdr = 0.005
cut_off_logFC = 2
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘None’，保存到change列
# 这里的change列用来设置火山图点的颜色
annotatedResults$Sig = ifelse(annotatedResults$FDR < cut_off_fdr & 
                   abs(annotatedResults$logFC) >= cut_off_logFC, 
                 ifelse(annotatedResults$logFC> cut_off_logFC ,'over-expressed','under-expressed'),'notDE')
table(annotatedResults$Sig)

library(ggplot2)

ggplot(annotatedResults, aes(x = logFC, y = -log10(FDR), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#d2dae2", "#ff4757","#546de5"))+
  #scale_color_manual(values=c("#ff4757", "#546de5","#d2dae2"))+
  # 辅助线
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_fdr),
             lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",
       y="-log10 (FDR)")+
  theme_bw()+
  ggtitle("DEG Volcano Plot")+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
