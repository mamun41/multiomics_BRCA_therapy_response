## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(123) # for reproducibility, remove for normal use

mydata <- read.csv("C:/Users/chenyuzhou/Desktop/feature.csv")
# Read the xlsx file

DEG <- read.csv("C:/Users/chenyuzhou/Desktop/diablo/DEG.csv")
microenvironment <- read.csv("C:/Users/chenyuzhou/Desktop/diablo/microenvironment.csv")
molecular_markers <- read.csv("C:/Users/chenyuzhou/Desktop/diablo/molecular_markers.csv")
somatic_mutations <- read.csv("C:/Users/chenyuzhou/Desktop/diablo/somatic_mutations.csv")

DEG <- scale(DEG)
molecular_markers <- scale(molecular_markers)
microenvironment <- scale(microenvironment)
somatic_mutations <- scale(somatic_mutations)

data = list(DEG=DEG,molecular_markers=molecular_markers,somatic_mutations=somatic_mutations,microenvironment=microenvironment)

lapply(data, dim) # check their dimensions

Y = as.factor(mydata[,c('RCB.category')]) # set the response variable as the Y dataframe
Y = as.factor(mydata[,c('resp.pCR')])
summary(Y)

pls1 <- spls(data[["DEG"]], data[["microenvironment"]]) 
pls2 <- spls(data[["DEG"]], data[["molecular_markers"]])
pls3 <- spls(data[["DEG"]], data[["somatic_mutations"]]) 
pls4 <- spls(data[["microenvironment"]], data[["molecular_markers"]])
pls5 <- spls(data[["microenvironment"]], data[["somatic_mutations"]]) 
pls6 <- spls(data[["molecular_markers"]], data[["somatic_mutations"]])

cor(pls1$variates$X, pls1$variates$Y) 
cor(pls2$variates$X, pls2$variates$Y) 
cor(pls3$variates$X, pls3$variates$Y)
cor(pls4$variates$X, pls4$variates$Y) 
cor(pls5$variates$X, pls5$variates$Y) 
cor(pls6$variates$X, pls6$variates$Y) 


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0
design[3,4] <- 0
design[4,3] <- 0
design

basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 10, design = design)


perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
perf.diablo$choice.ncomp$WeightedVote 

test.keepX = list (
  DEG=seq(8, 13, 1),
  molecular_markers=seq(9, 14, 1),
  immune_cell=seq(3,8,1),
  mutation=seq(9, 14, 1))

# run the feature selection tuning
tune.data = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")
saveRDS(tune.data, "selected_model.data")
library(reshape)
tune.data <- readRDS("selected_model.data")

names(tune.data[["choice.keepX"]]) [3]='somatic_mutations'
names(tune.data[["choice.keepX"]]) [4]='microenvironment'


list.keepX = tune.data$choice.keepX # set the optimal values of features to retain
list.keepX

final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)

selectVar(final.diablo.model, block = 'DEG', comp = 1)$DEG$name
selectVar(final.diablo.model, block = 'microenvironment', comp = 1)$immune_cell$name
selectVar(final.diablo.model, block = 'molecular_markers', comp = 1)$molecular_markers$name
selectVar(final.diablo.model, block = 'somatic_mutations', comp = 1)$mutation$name

plotDiablo(final.diablo.model, ncomp = 1)

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen','blue'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5,comp=1)

network(final.diablo.model, blocks = c(1,2,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen','blue'), cutoff = 0.7)

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15,14), cex = c(2,2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen','blue'))

cimDiablo(final.diablo.model,size.legend = 0.9,legend.position = "bottomright",)

perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate

auc.splsda = auroc(final.diablo.model, roc.block = "DEG", 
                   roc.comp = 1, print = FALSE)

auc.splsda = auroc(final.diablo.model, roc.block = "microenvironment", 
                   roc.comp = 1, print = FALSE)

auc.splsda = auroc(final.diablo.model, roc.block = "molecular_markers", 
                   roc.comp = 1, print = FALSE)

auc.splsda = auroc(final.diablo.model, roc.block = "somatic_mutations", 
                   roc.comp = 1, print = FALSE)

plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median')

