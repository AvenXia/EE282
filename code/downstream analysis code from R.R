library(Seurat)
library(dplyr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(monocle)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(Augur)
#Part1: Create Seurat Objects
plan("multiprocess", workers = 4)
plan()
dir <- c("data/GSE137478_RAW/control", "data/GSE137478_RAW/sample")                        
samples_name = c('GEO_control','GEO_sample')
scRNAlist_BQC <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist_BQC[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                           min.cells=3, min.features = 200)
  scRNAlist_BQC[[i]] <- RenameCells(scRNAlist_BQC[[i]], add.cell.id = samples_name[i])
  if(T){    
    scRNAlist_BQC[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist_BQC[[i]], pattern = "^mt-")
  }
  if(T){
    scRNAlist_BQC[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist_BQC[[i]], pattern = "^Rb[sl]")
  }
}

names(scRNAlist_BQC) <- samples_name
saveRDS(scRNAlist_BQC, file = "result_geo/scRNAlist_BQC.rds")
#Part2: Quality Control of sn-RNA data
load(file = "Rdata_result/p4-40/scRNAlist_BQC.Rdata")
scRNA_BQC <- merge(scRNAlist_BQC[[1]], scRNAlist_BQC[2:length(scRNAlist_BQC)])
png("result_geo/QC_before.png", width = 1800  , height = 600)
VlnPlot(scRNA_BQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
scRNA_AQC <- subset(scRNA_BQC,subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 10)
png("result_geo/QC_after.png", width = 1800  , height = 600)
VlnPlot(scRNA_AQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
scRNA_AQC <- scRNA_AQC[substr(rownames(scRNA_AQC),1,3) != "mt-"]
png("result_geo/QC_after_substr.png", width = 1800  , height = 600)
VlnPlot(scRNA_AQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
scRNA_AQC[["percent.mt"]] <- PercentageFeatureSet(scRNA_AQC, pattern = "^mt-")
png("result_geo/QC_after_substr_featureset.png", width = 1800  , height = 600)
VlnPlot(scRNA_AQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
save(scRNA_AQC, file = "Rdata_result/p4-40/scRNA_AQC.Rdata")
#Part3: Integration, PCA and clustering
scRNA_AQC <- readRDS(file = "data/scRNA_AQC_separateqc.rds")
scRNAlist <- SplitObject(scRNA_AQC, split.by = "orig.ident")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = scRNAlist)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features, reduction = "rpca", k.anchor = 20)
scrna_combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(scrna_combined) <- "integrated"
scrna_combined <- ScaleData(scrna_combined, verbose = FALSE)
scrna_combined <- RunPCA(scrna_combined, npcs = 30, verbose = FALSE)
scrna_combined <- RunUMAP(scrna_combined, reduction = "pca", dims = 1:30, umap.method = "umap-learn")
scrna_combined <- FindNeighbors(scrna_combined, reduction = "pca", dims = 1:30)
scrna_combined <- FindClusters(scrna_combined, resolution = 0.9)
pdf("result/umap.pdf", width = 12, height = 12)
DimPlot(scrna_combined, reduction = "umap", label = TRUE)
dev.off()
saveRDS(scrna_combined, file = "data/scrna_combined.rds")
#Part4: Identifying marker gene expression pattern by featureplots
scrna_combined <- readRDS(file = "data/scrna_combined.rds")
DefaultAssay(scrna_combined) <- "integrated"
pdf("result/UMAP/p4-40/filter.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/neuron.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Snap25", "Syt1", "Stmn2", "Rbfox3"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/glu.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Slc17a6", "Slc1a2", "Slc17a7"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Slc32a1", "Slc6a1","Gad1", "Gad2"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/astrocyte.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Aqp4", "Agt", "Sox9", "Gfap"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/endo.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c('Cldn5','Slc38a5', 'Pecam1', 'Icam2'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/opc.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c("Pdgfra", "Cspg4"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/ol.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Apod","Pmp22", "Trf", "Klk6"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/microglia.pdf", width = 12, height = 12) 
FeaturePlot(scrna_combined, features = c("C1qa", "Cx3cr1", "Spi1", "Tmem119"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/pericyte.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Pdgfrb"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/nol.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c("Tcf7l2", "Casr"), slot = "data") 
dev.off()
pdf("result/UMAP/p4-40/mol.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c("Mal", "Mog", "Plp1", "Serinc5"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/smc.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Acta2', 'Tagln', 'Myl9', 'Mylk'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/tan.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Nes', 'Lhx2', "Col23a1"), slot = "data")
dev.off()

pdf("result/UMAP/p4-40/gaba_pvalb.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Pvalb', 'Myo1e', 'Cemip', 'Nek7'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba_sst.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Sst', 'Unc13c', 'Lama4', 'Ptpru'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba_vip.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c('Vip', 'Sema5b'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba_lamp5.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Lamp5', 'Sv2c', 'Ndnf', 'Crispld1'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba_stac.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c('Stac', 'Cobll1'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/gaba_frem1.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Frem1', 'Megf10', 'March3'), slot = "data")
dev.off()

pdf("result/UMAP/p4-40/L2_3.pdf", width = 18, height = 18)
FeaturePlot(scrna_combined, features = c('Cux2', 'Ccbe1', 'Lhx2', 'Tle1', "Mdgal", "Otof"), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L4.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Rorb', 'Whrn', 'Rspo1', 'Kcnip2'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L5ET.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Fezf2'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L5IT.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Deptor', 'Rorb', 'Fezf2', 'Slc30a3'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L5NP.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Trhr', 'Slc17a8', 'Rapgef3', 'Sla2'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L5PT.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Bcl6'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L6CT.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c('Foxp2', 'Syt6'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L6IT.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Cdh9', 'Osr1', 'Sulf1', 'Slc30a3'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L6ITCAR3.pdf", width = 12, height = 12)
FeaturePlot(scrna_combined, features = c('Car3', 'Oprk1', 'Nr2f2'), slot = "data")
dev.off()
pdf("result/UMAP/p4-40/L6B.pdf", width = 12, height = 6)
FeaturePlot(scrna_combined, features = c('Nxph4', 'Cplx3'), slot = "data")
dev.off()
#Part5: 
celltype <- c("inhibitory neurons", "excitatory neurons", "oligodendrocytes", "inhibitory neurons", "oligodendrocytes", "excitatory neurons",
               "inhibitory neurons", "excitatory neurons", "endothelial cells", "excitatory neurons", "Oligodendrocyte precursor cells", "astrocytes", "excitatory neurons", 
               "inhibitory neurons", "excitatory neurons", "pericytes", "microglia", "excitatory neurons", "inhibitory neurons", "Vascular smooth muscle cells", "excitatory neurons", "excitatory neurons" )
Idents(scrna_combined) <- "seurat_clusters"
names(celltype) <- levels(scrna_combined)
scrna_combined <- RenameIdents(scrna_combined, celltype)
scrna_combined$celltype <- Idents(p40)
png("umap/umapcelltype.png", width = 1200  , height = 1200)
DimPlot(p40, reduction = "umap", pt.size = 0.5, label = TRUE)+ NoLegend()
dev.off()
##plotting dotplot of marker gene expression
scrna_combined$f2<-factor(scrna_combined@active.ident ,levels =c('Excitatory neurons', 'Inhibitory neurons',
                                                                 'Astrocytes','Endothelial cells','OPC',
                                                                 'VLMC',"nOL", 'OL','Microglia','Pericytes') )
Idents(scrna_combined) <- 'f2'
DefaultAssay(scrna_combined) <- 'RNA'
markers <- c('Stmn2',
             'Rbfox3',
             'Slc17a7',
             'Gad1',
             'Gad2',
             'Aqp4',
             'Gfap',
             'Cldn5',
             'Myh11',
             'Pdgfra',
             'Cspg4',
             'Tcf7l2',
             'Apod',
             'Trf',
             'C1qa',
             'Cx3cr1',
             'Pdgfrb')
pdf("result/dotplot/dotplot_celltype2.pdf", width = 17, height = 12)
DotPlot(scrna_combined, features = markers, group.by = "f2")+coord_flip()
dev.off()
#Part6: Re-clustering of ExN
scrna_combined_glu <- subset(scrna_combined_unfiler, ident=c(
  
))
DefaultAssay(scrna_combined_glu) <- 'RNA'
scRNAlist <- SplitObject(scrna_combined_glu, split.by = "orig.ident")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = scRNAlist)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features, reduction = "rpca", k.anchor = 20)
scrna_combined_glu <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(scrna_combined_glu) <- "integrated"
scrna_combined_glu <- ScaleData(scrna_combined_glu, verbose = FALSE)
scrna_combined_glu <- RunPCA(scrna_combined_glu, npcs = 15, verbose = FALSE)
scrna_combined_glu <- FindNeighbors(scrna_combined_glu, reduction = "pca", dims = 1:15)
scrna_combined_glu <- FindClusters(scrna_combined_glu, resolution = 0.6)
scrna_combined_glu <- RunUMAP(scrna_combined_glu, reduction = "pca", dims = 1:15)
pdf("result/UMAP/glu/umap.pdf", width = 15, height = 15)
DimPlot(scrna_combined_glu, reduction = "umap", label = TRUE)
dev.off()
DefaultAssay(scrna_combined_glu) <- 'RNA'
pdf("result/UMAP/glu/L2_3.pdf", width = 15, height = 18)
FeaturePlot(scrna_combined_glu, features = c('Cux2', 'Ccbe1', 'Lhx2', 'Tle1', "Mdga1", "Otof"), slot = "data")
dev.off()
celltype <- c('L2/3','L4','L2/3', 'L6CT', 'L4', 'L6IT', 'L6IT', 'L6CT', 'L5PT', 'L6IT', 'L5IT', 'L6b', 'C13', 'L4','L5NP','L4','L6CT','L6IT','C19','C20','L5PT','L5IT')
Idents(scrna_combined_glu) <- "seurat_clusters"
names(celltype) <- levels(scrna_combined_glu)
scrna_combined_glu <- RenameIdents(scrna_combined_glu, celltype)
scrna_combined_glu$celltype <- Idents(scrna_combined_glu)
Idents(scrna_combined) <- "celltype"
Idents(scrna_combined) <- "seurat_clusters"
pdf("result/UMAP/glu/umap_final_resol0.9_dim30_celltype.pdf", width = 14, height = 12)
DimPlot(scrna_combined_glu, reduction = "umap", label = TRUE, label.size = 6)
dev.off()
pdf("result/UMAP/glu/umap_final_resol0.9_dim30_celltype_split.pdf", width = 42, height = 12)
DimPlot(scrna_combined_glu, reduction = "umap", label = TRUE, label.size = 6,split.by = "orig.ident")
dev.off()
saveRDS(scrna_combined_glu, file = "result/scrna_combined_glu.rds")
augur <-  calculate_auc(scrna_combined_glu, cell_type_col = "celltype", label_col = "orig.ident")
#Part7: Differential gene expression analysis
GLU <- subset(scrna_combined, ident = "Excitatory neurons")
Idents(GLU) <- "orig.ident"
GLU_marker <- FindAllMarkers(GLU, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, only.pos = T) 
GLU_marker <- GLU_marker[,c(7,1:6)]
write.csv(GLU_marker,file = "result/DGE_Glu/GLU_allmarker.csv", row.names=F)
GLU_marker <- read.csv("result/DGE_Glu/csv/GLU_allmarker.csv")
GLU_marker_SE <- subset(GLU_marker, subset = cluster == "SE")
ego_GO <- list()
ego_GO_BP <- enrichGO(gene          = GLU_marker_SE$gene,
                      OrgDb         = 'org.Mm.eg.db',
                      keyType       = 'SYMBOL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
ego_GO_BP@result$Description <- substring(ego_GO_BP@result$Description,1,70) 
ego_GO_bp <- data.frame(ego_GO_BP)
write.csv(ego_GO_bp, file = "result/DGE_Glu/GO_BP.csv")
ego_GO_bp <- ego_GO_bp[order(ego_GO_bp$p.adjust),]
ego_GO_bp_top30 <- ego_GO_bp[1:30, ]
ego_GO_bp_top30 <- ego_GO_bp_top30[order(ego_GO_bp_top30$p.adjust,decreasing = T),]
ego_GO_bp_top30$Description <- factor(ego_GO_bp_top30$Description, levels = ego_GO_bp_top30$Description)
write.csv(ego_GO_bp_top30, file = "result/DGE_Glu/GO_BP_top30.csv")
pdf("result/DGE_Glu/GO_BP_top30.pdf", width = 12, height = 12)
ggplot(ego_GO_bp_top30, aes(y=Description, x = Count))+
  geom_bar(stat = "identity", width = 0.8, fill="salmon1")+
  theme_bw()+
  labs(x = "Num of genes", y = "Biological Process")+
  theme_classic(base_line_size = 0.6)+
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13), 
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 13,color="black"))
dev.off()

ego_GO_MF <- enrichGO(gene          = GLU_marker_SE$gene,
                      OrgDb         = 'org.Mm.eg.db',
                      keyType       = 'SYMBOL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
ego_GO_MF@result$Description <- substring(ego_GO_MF@result$Description,1,70) 
ego_GO_mf <- data.frame(ego_GO_MF)
write.csv(ego_GO_mf, file = "result/DGE_Glu/GO_MF.csv")
ego_GO_mf <- ego_GO_mf[order(ego_GO_mf$p.adjust),]
ego_GO_mf_top30 <- ego_GO_mf[1:30, ]
ego_GO_mf_top30 <- ego_GO_mf_top30[order(ego_GO_mf_top30$p.adjust,decreasing = T),]
ego_GO_mf_top30$Description <- factor(ego_GO_mf_top30$Description, levels = ego_GO_mf_top30$Description)
write.csv(ego_GO_mf_top30, file = "result/DGE_Glu/GO_MF_top30.csv")
pdf("result/DGE_Glu/GO_MF_top30.pdf", width = 12, height = 12)
ggplot(ego_GO_mf_top30, aes(y=Description, x = Count))+
  geom_bar(stat = "identity", width = 0.8, fill="salmon1")+
  theme_bw()+
  labs(x = "Num of genes", y = "Molecular Function")+
  theme_classic(base_line_size = 0.6)+
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13), 
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 13,color="black"))
dev.off()

Idents(scrna_combined_glu) <- "orig.ident"
top10 <- read.csv("result/DGE_Glu/csv/GLU_allmarker_heat.csv")
matric <- AverageExpression(scrna_combined_glu, assays = "RNA")[[1]]
matric <- matric[rownames(matric)%in%top10$gene,]
choose_matrix = t(scale(t(log2(matric+1))))
pdf("result/DGE_Glu/heatmap_genelist.pdf", width = 12, height = 12)
pheatmap(choose_matrix,  scale = "none",  clustering_method = "mcquitty", border_color = NA, cluster_cols = F, cluster_rows = T, fontsize_row = 10, cellheight = 8, cellwidth = 75)
dev.off()

#Part7: Pseudotime analysis of Layer II/III ExN
GLU <- SplitObject(GLU, split.by = 'celltype')
mycds <- list()
deg <- list()
name2 <- c("L4 ordering_gene.pdf", "L6CT ordering_gene.pdf", "L2_3 ordering_gene.pdf", "L6IT ordering_gene.pdf", "L5PT ordering_gene.pdf", "L5NP ordering_gene.pdf", "L6b ordering_gene.pdf", "L5IT ordering_gene.pdf")
name1 <- c("L4 deg.csv", "L6CT deg.csv", "L2_3 deg.csv", "L6IT deg.csv", "L5PT deg.csv", "L5NP deg.csv", "L6b deg.csv", "L5IT deg.csv")
name3 <- c("L4 state_tra.pdf", "L6CT state_tra.pdf", "L2_3 state_tra.pdf", "L6IT state_tra.pdf", "L5PT state_tra.pdf", "L5NP state_tra.pdf", "L6b state_tra.pdf", "L5IT state_tra.pdf")
name4 <- c( "L4 orig_tra.pdf", "L6CT orig_tra.pdf", "L2_3 orig_tra.pdf", "L6IT orig_tra.pdf", "L5PT orig_tra.pdf", "L5NP orig_tra.pdf", "L6b orig_tra.pdf", "L5IT orig_tra.pdf")
name5 <- c("L4 pseu_tra.pdf", "L6CT pseu_tra.pdf", "L2_3 pseu_tra.pdf", "L6IT pseu_tra.pdf", "L5PT pseu_tra.pdf", "L5NP pseu_tra.pdf", "L6b pseu_tra.pdf", "L5IT pseu_tra.pdf")
name6 <- c( "L4 mycds.rds", "L6CT mycds.rds", "L2_3 mycds.rds", "L6IT mycds.rds", "L5PT mycds.rds", "L5NP mycds.rds", "L6b mycds.rds",  "L5IT mycds.rds")
for (i in 1:length(GLU)){
  Idents(GLU[[i]]) <- 'orig.ident'
  expr_data <- as(as.matrix(GLU[[i]]@assays$RNA@counts), 'sparseMatrix')
  p_data <- GLU[[i]]@meta.data
  p_data$orig.ident <- GLU[[i]]@active.ident
  pd <- new('AnnotatedDataFrame', data = p_data)
  f_Data <- data.frame(gene_short_name = row.names(GLU[[i]]), row.names = row.names(GLU[[i]]))
  fd <- new('AnnotatedDataFrame', data = f_Data)
  mycds[[i]] <- newCellDataSet(expr_data,
                               phenoData = pd,
                               featureData = fd,
                               expressionFamily = negbinomial.size())
  mycds[[i]] <- estimateSizeFactors(mycds[[i]])
  mycds[[i]] <- estimateDispersions(mycds[[i]])
  mycds[[i]] <- detectGenes(mycds[[i]], min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(mycds[[i]]),
                                      num_cells_expressed >= 10))
  diff_test_res <- differentialGeneTest(mycds[[i]][expressed_genes,],
                                        fullModelFormulaStr = "~orig.ident")
  DEG <- subset(diff_test_res, qval < 0.01)
  deg[[i]] <- DEG[order(DEG$qval,decreasing = F),]
  write.csv(deg[[i]], file = name1[[i]], col.names = T, row.names = F, sep = "\t", quote = F)
  ordering_genes <- row.names (deg[[i]])
  mycds[[i]] <- setOrderingFilter(mycds[[i]], ordering_genes)
  pdf(name2[[i]], width = 12, height = 12)
  plot_ordering_genes(mycds[[i]])
  dev.off()
  mycds[[i]] <- reduceDimension(mycds[[i]], max_components = 2,
                                method = 'DDRTree')
  mycds[[i]] <- orderCells(mycds[[i]])
  pdf(name3[[i]], width = 12, height = 12)
  plot_cell_trajectory(mycds[[i]], color_by = 'State', size = 1, show_backbone = T)
  dev.off()
  pdf(name4[[i]], width = 12, height = 12)
  plot_cell_trajectory(mycds[[i]], color_by = 'orig.ident', size = 1, show_backbone = T)
  dev.off()
  pdf(name5[[i]], width = 12, height = 12)
  plot_cell_trajectory(mycds[[i]], color_by = 'Pseudotime', size = 1, show_backbone = T)
  dev.off()
  saveRDS(mycds[[i]], file = name6[[i]])
}

mycds_L2 <- readRDS("result/PSEUDO/L2_3 mycds.rds")
mycds_L4 <- readRDS("result/PSEUDO/L4 mycds.rds")
mycds_L5IT <- readRDS("result/PSEUDO/L5IT mycds.rds")
mycds_L5NP <- readRDS("result/PSEUDO/L5NP mycds.rds")
mycds_L5PT <- readRDS("result/PSEUDO/L5PT mycds.rds")
mycds_L6b <- readRDS("result/PSEUDO/L6b mycds.rds")
mycds_L6IT <- readRDS("result/PSEUDO/L6IT mycds.rds")
mycds_L6CT <- readRDS("result/PSEUDO/L6CT mycds.rds")

BEAM_res_L2 <- BEAM(mycds_L2, branch_point = 2, cores = 8)
BEAM_res_L4 <- BEAM(mycds_L4, branch_point = 1, cores = 8)
BEAM_res_L5NP <- BEAM(mycds_L5NP, branch_point = 2, cores = 8)
BEAM_res_L6IT <- BEAM(mycds_L6IT, branch_point = 1, cores = 8)

saveRDS(BEAM_res_L2, file = 'sig2.rds')
saveRDS(BEAM_res_L4, file = 'sig4.rds')
saveRDS(BEAM_res_L5NP, file = 'sig5.rds')
saveRDS(BEAM_res_L6IT, file = 'sig6.rds')

#Part8: Reclustering analysis of microglia
scrna_combined_micro <- subset(scrna_combined_unfiler, ident=c(
  '4','41'
))
DefaultAssay(scrna_combined_micro) <- 'RNA'
scRNAlist <- SplitObject(scrna_combined_micro, split.by = "orig.ident")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = scRNAlist)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features, reduction = "rpca", k.anchor = 20)
scrna_combined_micro <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(scrna_combined_micro) <- "integrated"
scrna_combined_micro <- ScaleData(scrna_combined_micro, verbose = FALSE)
scrna_combined_micro <- RunPCA(scrna_combined_micro, npcs = 15, verbose = FALSE)
scrna_combined_micro <- FindNeighbors(scrna_combined_micro, reduction = "pca", dims = 1:15)
scrna_combined_micro <- FindClusters(scrna_combined_micro, resolution = 0.6)
scrna_combined_micro <- RunUMAP(scrna_combined_micro, reduction = "pca", dims = 1:15)
pdf("result/UMAP/micro/umap_micro.pdf", width =8, height = 8)
DimPlot(scrna_combined_micro, reduction = "umap", label = TRUE)
dev.off()
scrna_combined_micro <- subset(scrna_combined_micro, ident = c("0","1","2","3","4","5"))
saveRDS(scrna_combined_micro, file = "result/scrna_combined_micro.rds")
DefaultAssay(scrna_combined_micro) <- 'RNA'
celltype <- c("AMG1", "Homeostatic microglia", "AMG2", "AMG3", "AMG4", "Macrophage")
Idents(scrna_combined_micro) <- "seurat_clusters"
names(celltype) <- levels(scrna_combined_micro)
scrna_combined_micro <- RenameIdents(scrna_combined_micro, celltype)
scrna_combined_micro$celltype <- Idents(scrna_combined_micro)
Idents(scrna_combined_micro) <- "celltype"
Idents(scrna_combined_micro) <- "seurat_clusters"
pdf("result/UMAP/micro/micro_final_resol0.9_dim30_celltype.pdf", width = 10, height = 8)
DimPlot(scrna_combined_micro, reduction = "umap", label = TRUE, label.size = 5)
dev.off()
pdf("result/UMAP/micro/split_micro.pdf", width = 30, height = 10)
DimPlot(scrna_combined_micro, reduction = "umap", split.by = "orig.ident", label = T, label.size = 8)
dev.off()

DefaultAssay(scrna_combined_micro) <- "RNA"
marker_micro <- FindAllMarkers(scrna_combined_micro, slot = "data", logfc.threshold = 0.25, min.pct = 0.1, only.pos = T)
marker_micro <- marker_micro[,c(7,1:6)]
write.csv(marker_micro, "result/cluster/micro/marker_micro.csv", row.names = F)
marker_micro <- read.csv(file = "result/cluster/micro/marker_micro.csv")
AMG3 <- subset(marker_micro, cluster == "3")

AMG3 <- AMG3[order(AMG3$difference,decreasing = T),]
c13_top100 <- AMG3[1:200,]
GO_c13_BP <- enrichGO(gene          = c13_top100$gene,
                      OrgDb         = 'org.Mm.eg.db',
                      keyType       = 'SYMBOL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
GO_c13_BP@result$Description <- substring(GO_c13_BP@result$Description,1,70) 
GO_c13_bp <- data.frame(GO_c13_BP)
write.csv(GO_c13_bp, file = "result/cluster/micro/GO_AMG3_BP.csv")
GO_c13_bp <- GO_c13_bp[order(GO_c13_bp$p.adjust),]
GO_c13_bp_top30 <- GO_c13_bp[1:30, ]
GO_c13_bp_top30 <- GO_c13_bp_top30[order(GO_c13_bp_top30$p.adjust,decreasing = T),]
GO_c13_bp_top30$Description <- factor(GO_c13_bp_top30$Description, levels = GO_c13_bp_top30$Description)
write.csv(GO_c13_bp_top30, file = "result/cluster/micro/GO_AMG3_bp_top30.csv")
pdf("result/cluster/micro/GO_AMG3_bp_top30.pdf", width = 12, height = 12)
ggplot(GO_c13_bp_top30, aes(y=Description, x = Count))+
  geom_bar(stat = "identity", width = 0.8, fill="salmon1")+
  theme_bw()+
  labs(x = "Num of genes", y = "Biological Process")+
  theme_classic(base_line_size = 0.6)+
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13), 
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 13,color="black"))
dev.off()
GO_c13_MF <- enrichGO(gene          = AMG3_top200$gene,
                      OrgDb         = 'org.Mm.eg.db',
                      keyType       = 'SYMBOL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
GO_c13_MF@result$Description <- substring(GO_c13_MF@result$Description,1,70) 
GO_c13_mf <- data.frame(GO_c13_MF)
write.csv(GO_c13_mf, file = "result/cluster/micro/GO_AMG3_mf.csv")
GO_c13_mf <- GO_c13_mf[order(GO_c13_mf$p.adjust),]
GO_c13_mf_top30 <- GO_c13_mf[1:30, ]
GO_c13_mf_top30 <- GO_c13_mf_top30[order(GO_c13_mf_top30$p.adjust,decreasing = T),]
GO_c13_mf_top30$Description <- factor(GO_c13_mf_top30$Description, levels = GO_c13_mf_top30$Description)
write.csv(GO_c13_mf_top30, file = "result/cluster/micro/GO_AMG3_mf_top30.csv")
pdf("result/cluster/micro/GO_AMG3_mf_top30.pdf", width = 12, height = 12)
ggplot(GO_c13_mf_top30, aes(y=Description, x = Count))+
  geom_bar(stat = "identity", width = 0.8, fill="salmon1")+
  theme_bw()+
  labs(x = "Num of genes", y = "Molecular Function")+
  theme_classic(base_line_size = 0.6)+
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13), 
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 13,color="black"))
dev.off()
proliferation <- c("Pola1", "Mki67", "Mcm10", "Mcm6", "Mcm9")
plot <- list()
for (i in 1:length(proliferation)){
  plot[[i]] <-  FeaturePlot(scrna_combined_micro, features = proliferation[i]) 
}
pdf("result/cluster/micro/AGM3_MARKER.pdf", width = 12, height = 12)
for (i in 1:length(proliferation)){
  print(plot[[i]])
}
dev.off()


