library(dplyr)
library(Seurat)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(ggplot2)

setwd("D:\\AIDA_project\\data\\H136_result")
human.data <- Read10X(data.dir = "../H136_filtered_feature_bc_matrix/") 
human_obj <- CreateSeuratObject(counts = human.data, project = "human_singlecell", min.cells = 3, min.features = 200)  
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
human_obj[["percent.mt"]] <- PercentageFeatureSet(human_obj, pattern = "^MT-")
# Visualize QC metrics as a violin plot
pdf("01.vlnplot_QC.pdf",width=6,height = 4)
VlnPlot(human_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02.feature_plot.pdf",width=12,height = 4)
plot1 <- FeatureScatter(human_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

human_obj <- subset(human_obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalizing the data
human_obj <- NormalizeData(human_obj, normalization.method = "LogNormalize", scale.factor = 10000)
human_obj <- NormalizeData(human_obj)

# Identification of highly variable features (feature selection)
human_obj <- FindVariableFeatures(human_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(human_obj), 10)

# plot variable features with and without labels
pdf("03.top10_variable_genes.pdf",width=12,height = 4)
plot1 <- VariableFeaturePlot(human_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Scaling the data
all.genes <- rownames(human_obj)
human_obj <- ScaleData(human_obj, features = all.genes)

# Perform linear dimensional reduction
human_obj <- RunPCA(human_obj, features = VariableFeatures(object = human_obj))
# Examine and visualize PCA results a few different ways
print(human_obj[["pca"]], dims = 1:5, nfeatures = 5)
#DimPlot(human_obj, reduction = "pca")
#VizDimLoadings(human_obj, dims = 1:2, reduction = "pca")

# Determine the ¡®dimensionality¡¯ of the dataset
ElbowPlot(human_obj)

# Cluster the cells
human_obj <- FindNeighbors(human_obj, dims = 1:20)
human_obj <- FindClusters(human_obj, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP/tSNE)
reticulate::py_install(packages ='umap-learn')
human_obj <- RunUMAP(human_obj, dims = 1:20)
human_obj <- RunTSNE(human_obj, dims = 1:20)
saveRDS(human_obj,sep = "")
#saveRDS(human_obj, file = paste(parameter_workdir,"/","after_filter.rds",sep = ""))

# individual clusters
pdf("04.UMAP_plot.pdf",width = 6,height = 4)
DimPlot(human_obj, reduction = "umap",)
dev.off()

pdf("05.TSNE_plot.pdf",width = 6,height = 4)
DimPlot(human_obj, reduction = "tsne")
dev.off()

##### annotation ######
# find markers for every cluster compared to all remaining cells, report only the positive
pbmc.markers <- FindAllMarkers(human_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
pdf("06.heatmap.pdf",width = 8,height = 16)
DoHeatmap(human_obj, features = top10$gene) + NoLegend()
dev.off()

## DEG markers by different clusters

if(!dir.exists("DEGcluster_markers")){dir.create("DEGcluster_markers")}
setwd("DEGcluster_markers")
GeneBoxplot <- function(Object = Object, Object_Mat = Object_Mat, Feature = Feature){
  require(ggplot2)
  TargetGeneMat <- Object_Mat[Feature,]
  TmpMat <- data.frame(Cluster = factor(x = as.character(Idents(Object)), levels = levels(Idents(Object))), Expression = Object_Mat[Feature,])
  p <- ggplot(TmpMat, aes(x = Cluster, y = Expression)) +
    geom_boxplot(show.legend = F, outlier.size = 1, outlier.fill = NULL , aes(fill= Cluster)) +
    labs(y = Feature) +
    theme_bw() +
    theme(axis.text.y = element_text(colour = "black"),
          panel.grid=element_blank(),
          axis.text.x = element_text(angle = 90, colour = "black"),
          panel.border=element_blank(),
          axis.line=element_line(size=1,colour="black")) +
    scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Set1"))(clu_num))
  return(p)
}

seurat_data_mat <- as.matrix(human_obj@assays$RNA@data)
marker<-data.frame(Cluster=top20$cluster,Marker=top20$gene)
clu_num<-length(levels(human_obj@meta.data$seurat_clusters))

for(i in unique(marker$Cluster)){
    gene<-intersect(rownames(human_obj@assays$RNA@data),marker[which(marker$Cluster==i),"Marker"])
    if(!dir.exists(i) & length(gene)!=0 ){
        dir.create(i)
        setwd(i)
        for (j in gene) {
            p1<-GeneBoxplot(Object = human_obj, Object_Mat = seurat_data_mat, Feature = j)
            p2<-FeaturePlot(human_obj, features = j)+scale_fill_gradientn(colors = c("#969696","yellow","red"))
            p3<-VlnPlot(human_obj, features = j , pt.size = 0)+scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Set1"))(clu_num))+labs(x="Cluster")
            pdf(paste(i,"_",j,".pdf",sep = ""),width=14)
            p<-(p1/p3)-p2 + plot_layout(widths = c(1,2))
            print(p)
            dev.off()
        }
        setwd("../")
    }
}
setwd("../")

## GO/KEGG pathway 
library(ggplot2)
library(reshape2)
library(ggthemes)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Hs.eg.db)
ClusterMarker <-pbmc.markers
setwd("..")
if(!dir.exists("GOKEGG")){dir.create("04_GOKEGG")}
setwd("GOKEGG")

for (i in sort(unique(ClusterMarker$cluster))){
    mygenes<-as.character(as.character(ClusterMarker[which(ClusterMarker$cluster == i),"gene"]))
    #if(length(mygenes) > 10){
    mygenes.df <- bitr(mygenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    head(mygenes.df,2)
    ego <- enrichGO(gene = mygenes.df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "All",
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
    ego.df <- data.frame(ego)
    ego.S <- ego.df[which(ego.df$p.adjust <= 0.05), ]
    write.table(ego.S, paste(paste("Sign_GO-C",i,sep = "_"),"csv",sep = "."), sep = ",", row.names = F)
  
    pathway <- as.data.frame(ego)
    pathway = pathway[1:20,]
    pathway = pathway[order(pathway$pvalue, decreasing = T),]
    tt <- factor(pathway$Description, levels = unique(pathway$Description))
    pp = ggplot(pathway, aes(-1*log10(p.adjust), tt))
    pbubble = pp + geom_point(aes(size = Count, color = -1*log10(p.adjust)))+
        scale_colour_gradient(low = "blue",high = "red") +
        theme_bw() +scale_size(range = c(2, 10)) +
        labs(title = paste("GO-enrich_cluster",i,sep = ""))+
        ylab("") +
        xlim(min(-1*log10(pathway$p.adjust))*0.9, max(-1*log10(pathway$p.adjust))*1.1)+
        theme(axis.text.y = element_text(size = 20, family = "Helvetica", color = "black", angle = 0),
              axis.text.x = element_text(size = 15, family = "Helvetica", color = "black", angle = 0))+
        theme(legend.title = element_text(color = "black", size = 20, family = "Helvetica"))+
        theme(legend.text = element_text(color="azure4", size = 15, family = "Helvetica"))+
        theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
        theme(legend.key.size = unit(1.1,'cm'))
    ggsave(file = paste(paste("GO-C",i,sep="_"),"pdf",sep="."), plot = pbubble, width = 20, height = 10)
  #}
}

## Marker genes for every cell type
#marker<-read.csv("../pbmc_marker.csv")
# deal with csv (from vector to data.frame)
#data <- as.matrix(strsplit(marker$B_pan,split = ",")[[1]])
#data <- data.frame(id=c(1:nrow(data)),data)
#names(data)[2] <-name[1]
#name <- names(marker)
#for (i in 2:length(name)) {
 # a <- as.matrix(strsplit(as.character(marker[i]),split = ",")[[1]])
 # a <- data.frame(id=c(1:nrow(a)),a)
 # names(a)[2] <-name[i]
#  data <- merge(data,a,all = T,by="id")
#}
#write.csv(data,"markergene.csv")
library(tidyverse)
data<-read.csv("markergene.csv",na.strings = F)
data<-data[,-1]
if(!dir.exists("cell_type_marker")){dir.create("cell_type_marker")}
setwd("cell_type_marker")
for (i in 1:length(colnames(data))) {
  c<-unique(c(na.omit(data[,i])))
  DotPlot(human_obj,features = c)+RotatedAxis()+
    scale_x_discrete("")+scale_y_discrete("")
  ggsave(paste0(colnames(data)[i],"_marker_expression.pdf"),width = 8,height = 5)
}
