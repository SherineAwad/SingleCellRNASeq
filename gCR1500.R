library(Seurat)
library(dplyr)
library(Matrix)




mydata.data <- Read10X(data.dir = "Gout/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
mydata <- CreateSeuratObject(raw.data = mydata.data, project = "Cell Ranger/Seurat", min.cells = 3, min.features = 200)




mito.genes <- grep(pattern = "^MT-", x = rownames(x = mydata@data), value = TRUE)
length(mito.genes)

# check out the meta data
head(mydata@meta.data)

percent.mito <- Matrix::colSums(mydata@raw.data[mito.genes, ]) / Matrix::colSums(mydata@raw.data)

# add some more meta data
mydata <- AddMetaData(object = mydata,
                    metadata = percent.mito,
                    col.name = "percent.mito")

head(mydata@meta.data)

# plot number of genes, UMIs, and % mitochondria
png(filename ="gCR.mito.png",width = 600, height = 600, units = "px")

VlnPlot(object = mydata,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)


dev.off()


# check how many genes have at least one transcript in each cell
png("gCR.geneswith1tr.png")
at_least_one <- apply(mydata.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
dev.off()
png("gCR.sumexpression.png")
hist(colSums(mydata.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
dev.off()



mydata <- FilterCells(object = mydata,
                    subset.names = c("nGene"),
                    low.thresholds = c(1500),
                    high.thresholds = c(3500))

png("gCR1500.before.hist.png")
hist(colSums(mydata@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()

mydata <- NormalizeData(object = mydata,
                      normalization.method = "LogNormalize",
                      scale.factor = 1e4)


png("gCR1500.after.hist.png")
hist(colSums(mydata@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")

dev.off() 

mydata <- FindVariableGenes(object = mydata,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)
mydata <- ScaleData(object = mydata,vars.to.regress = c("nUMI") )

mydata <- RunPCA(object = mydata,pc.genes = mydata@var.genes,do.print = TRUE,pcs.print = 1:5,genes.print = 5)

png("gCR1500.pca.png")

PCAPlot(object = mydata, dim.1 = 1, dim.2 = 2)

dev.off()


png("gCR1500.dispersion.png")
mydata  <- FindVariableGenes(object = mydata,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)
dev.off ()

png("gCR1500.heatmap1.png") 

PCHeatmap(object = mydata,
          pc.use = 1,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE,size.x.use =6, size.y.use =6)
dev.off() 

png("gCR1500.heatmap2.png") 
PCHeatmap(object = mydata,
          pc.use = 1:12,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE,size.x.use =6, size.y.use =6)
dev.off() 

mydata <- FindClusters(object = mydata,
                     reduction.type = "pca",
                     dims.use = 1:10,
                     resolution = 0.6,
                     print.output = 0,
                     save.SNN = TRUE)
 
PrintFindClustersParams(object = mydata)


png("gCR1500.tsne.cluster.png")

mydata <- RunTSNE(object = mydata,
                dims.use = 1:10,
                do.fast = TRUE)
 
TSNEPlot(object = mydata, do.label = TRUE)

dev.off()

head(PCTopCells(object = mydata, pc.use = 1, num.cells = NULL, do.balanced = FALSE))
head(PCTopGenes(object =mydata, pc.use = 1, num.genes = 30, use.full = FALSE,do.balanced = FALSE) )

png("gCR1500.featureplot.pca.png")
FeaturePlot(object = mydata,
            features.plot = c("STMN2","RTN1","GAP43","DCX","INA","MAP1B"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

dev.off()


png("gCR1500.vnplot.pca.png") 
VlnPlot(object = mydata,
        features.plot = c("STMN2","RTN1","GAP43","DCX","INA","MAP1B"), 
        use.raw = TRUE,
        y.log = TRUE)

dev.off() 


#Find all markers
head(mydata.markers <- FindAllMarkers(object = mydata,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25))

png("gCR1500.featureplot.marker.png")
FeaturePlot(object = mydata,  
            features.plot = c("SPARCL1","SLC1A3","GPC3","VIM","EEF1A1","MGST1"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

dev.off()

png("gCR1500.vnplot.marker.png")
VlnPlot(object = mydata,
        features.plot = c("SPARCL1","SLC1A3","GPC3","VIM","EEF1A1","MGST1"), 
        use.raw = TRUE,
        y.log = TRUE)

dev.off()



 
