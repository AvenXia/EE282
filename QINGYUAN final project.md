---
title: Transcriptome-based approach to dissect cell-type specific molecular mechanisms for mice focal cortical dysplasia
Name: Qingyuan Xia
word counts: 2000
---

# Introduction
Focal cortical dysplasia (FCD) is a developmental abnormality of cortical architecture. It is characterized by disorganization of cerebral cortex and cortical dysgenesis (Marsan, E. and Baulac, S., 2018) . In this study, we sought to uncover the cell type-specific molecular mechanisms in the FCD-induced epilepsy development. Single-nucleus RNA sequencing of mice neocortex at middle (GS) and severe (SE) stages of epilepsy was performed. We screened the upregulated transcriptional factors that are highly expressed at the severe stage of the epilepsy. Additionally, the proliferation potentials of epilepsy brain were identified, where the microglia is one candidates contributing to the proliferation. Together these findings shed light on the alteration of transcriptomic levels from neuron and microglia. The abnormal neuron states were marked with upregulation of neurodevelopment-related transcriptional factors. 

# Materials and Methods
## Raw data to be used
The datasets to be utilized is the single nuclei RNA sequencing data of neocortex from three transgenetic mice from my lab. Three groups of sn-RNA datasets were collected at different stages of of epilepsy: wildtype (negative control),mild (GS stage) and severe (SE stage). 
## Upstream analysis by bash coding
To do the sequence assembly, annotation and gene counting of the raw data, the pipeline called Cell Ranger from 10X Genomics is applied. It can  perform sample demultiplexing, barcode processing, single cell 3' and 5' gene counting, V(D)J transcript sequence assembly and annotation, and Feature Barcode analysis from single cell data. Raw reads have been first mapped to the mm10 mouse reference genome and demultiplexed to generate a per-cell count matrix using CellRanger pipeline. By running cellranger count, we picked up "filtered_feature_bc_matrix" from the output dataset, which will be applied for downstream R analysis.
## Downstream analysis by R
The following downstream analysis is achieved by Seurat from R. The input data will be filtered_feature_bc_matrix by cell ranger. Cells of potentially low quality (unique feature counts > 6000 or total UMI < 700 or total UMI > 30,000, or percentage of mitochondrial genes > 10%) were removed from downstream analysis. 
### Identifying diverse types of cell populations across multiple datasets
Identifying the cell populations across multiple datasets from different mice brings about the unique problem because of experimental variations (batch effects and disease conditions). To minimize this, the reciprocal PCA method in Seurat package is utilized that projects each dataset into the others PCA space and constrains the anchors by the same mutual neighborhood requirement. Briefly, gene counts are normalized and 2000 highly variable genes are selected in each dataset based on the mean/variance calculations. The features that are repeatedly variable across datasets is selected for integration PCA. After each dataset is scaled and performed with PCA, anchors between individual data are identified for data integration. The integrated data is further scaled and centered followed by regular PCA on purpose of dimension reduction. PC1 to PC30 are chosen to construct a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space, and the edge weights between any two cells are modified based on Jaccard similarity. The cell clustering in the KNN graph is achieved by applying Louvain algorithm clustering (resolution = around 1). For microglia/neuron re-clustering, the clusters were extracted and then clustered based on the identical integrative clustering and DGE method mentioned above.
### Finding the highly expressed genes in different stages of epilepsy
To identify the highly expressed genes in different stages of epilepsy compared with wildtype mice model, clusters of neuron were extracted based on the canonical markers. The expression value of each gene in each epilepsy stage were compared against the wildtype mouse data using Wilcoxon rank sum test. Upregulated genes recommended with p_value < 0.05, log fold change > 0.25 and minimized pct_difference > 0.1 can be selected. The DEG genelists can be visualized with heatmap using packages pheatmap. For the upregulated genes screened via DGE analysis, gene ontology database (GO Biological Process, GO Molecular Function, and GO Cellular Component) at mouse genome informatics was referenced via R package enrichProfiler to find out enriched ontology (adjusted p_value < 0.1) from the DGE list. 
### Simulating temporal trajectories of gene expression patterns of specific cell type during epilepsy development
Pseudo-time analysis is conducted via R package Monocle to find out the linear gene expression alterations. First a specific data of one cell type cluster should be extracted. For monocle analysis, default parameters are used. Cells in WT and every stages of the seizures are ordered into pseudo-time trajectory using the top 2000 high variable genes via orderCells function. Ordering genes with qval < 0.1 for pseudo-time ordering are considered using the differentialGeneTest() function. Ordering genes are clustered based on the expression trajectory using the ward.D2 method from the R package pheatmap. The expression level tracking of ordering genes over pseudo-time can also be visualized by function plot_genes_in_pseudotime from monocle.

# Results
## Quality control of data
From the feature, count data of RNA, it was found that some data has a higher number of count and feature that are not close to the main cluster. These data have a high risk of being double droplet data. To remove low quality data, threholds of feature (6000) and count (30000) of RNA are set. After quality control, all data points are clustered to a similar range (Figure 1). 

## Description of the data
An object of class Seurat 
27422 features across 18672 samples within 2 assays 
Active assay: RNA (25422 features, 0 variable features)
 1 other assay present: integrated
 2 dimensional reductions calculated: pca, umap

## Creating a single cell atlas of mice neocortex at different stages of epilepsy
After the removal of low-quality nuclei and dataset integration, a total dataset of 18672 nuclei was obtained. By the method of integration, data from three samples can be distributed in a big map. Based on the well-established reference markers from previous single cell sequencing analysis findings, 9 cell types were clustered including glutamatergic neurons, GABAergic neurons, endothelial cells, vascular smooth muscle cells, microglia, astrocytes, oligodendrocytes, and oligodendrocyte precursor cells (OPCs) (Figure 2). 

The excitatory neuron can be further re-classifeid into subtypes based on different layers of neocortex. Thus, ExN is further extracted and re-clustered into layer 2 to 6. The excitatory neurons at different cortical layers participate in brain activities from different neuronal circuits. The distinct neuronal inputs from the circuits may exhibit diverse neuronal sensitivities to the epilepsy transition. Therefore, we planned to determine which cortical layers responded most significantly to the epilepsy. Noticeably, statistical differences in gene expression can largely be driven by the size of the cluster, with larger clusters having the ability to resolve more differentially expressed genes. Thus, a statistical method Augur was applied that can rank responsiveness of cell types in single cell data without bias from cluster size. It was found that the average cell type prioritization score across all cell types increased during epilepsy development, suggesting distinct transcriptomic transition at both GS and SE stages. We found that ExN on cortical layer II/III were the most perturbation-responsive cell type at both epilepsy stages (0.84 and 0.95 average AUC at GS and SE respectively) (Figure 3).

## Epilepsy stage specific upregulation of gene expression profiles in glutamatergic neurons
Glutamatergic neurons, as the major excitatory neurons, had upregulated mTOR signaling due to Pten knockout, which is also the major cell type contributing to abnormal brain activities. Thus, we focused on uncovering its alteration of transcriptome during different stages of epilepsy.

After the ExN cluster was extracted from the cell atlas, we find out highly expressed genes by comparing SE/GS group with control group, which is shown from table1 (SE VS CTR) and table2(GS VS CTR). From the analysis, we can list highly expressed genes of each group. The index pct.1 and pct.2 refer the proportion of cells expressing this gene. By calculating pct.1-pct.2, we can tell which genes have a higher specifically gene expression pattern. To find out the the enrichment of DEGs on which biological processes, we performed Gene ontology analysis. At GS stages, the top 100 highly expressed genes were enriched in the biological processes associated with synapse development. Yet diverse types of activity dependent genes (ADG) (e.g. Fos, Jun, Junb, Egr1 and Nr4a1) were upregulated at the SE stages, which were enriched in neuronal apoptosis and post-translational protein folding. Thus, ExN at GS stage indicates an overactivation of synapse organization while SE stage exhibits a distinct protein synthesis stress and a high risk of inducing irreversible neuron damage (Figure 4).

Transcriptional factors (TFs) are critical DNA-binding proteins regulating gene expression, which is largely involved in the reorganization of neuronal activities. It was hypothesized that upregulation of a series of SE-specific TFs promotes the severity of epilepsy. we selected genes that satisfied three criteria across all ExN: (1) Their expression at SE should be higher than that at GS stage; (2) the expression at epilepsy stage should not be lower than wildtype. As a result, a total set of 75 TFs was screened (Figure 4). 

We wanted to see whether the TFs screened above also revealed a continuous alteration of expression levels during epilepsy development. ExN at cortical layer II/III as the most sensitive neuronal types to condition alteration, so we are trying to find out the molecule dynamics in this cell type. We carried out pseudotime analysis to order cells by differentiation state, where wildtype occurred early while SE ends in pseudotime. Among all genes across ExN, 1923 had significantly different trajectories at SE compared with wildtype. We then clustered significantly different gene trajectories to identify those with similar expression patterns. We find a Cluster that was characterized by a continuous increased expression following pseudotime. We noted that from this cluster, Id2, Rai14 and Klf5 were presented (Figure 5). 

## An increased proliferative potential from microglia
The activation of microglia is reported to be highly related to neuropathology. To further study the cell type specific mechanism of microglia, we performed the re-clustering of microglia to find out proliferative sub-population. Six subpopulations of microglia were observed, including one homeostasis microglia, four activated microglia (ACM1-4) and M2 macrophage. The homeostasis microglia were defined by P2ry12 and Trem2, and mainly exist at GS stage. The proportion of ACM1 and ACM2 is higher at SE and were defined by activity-dependent genes such as Hspa1a and Hspa1b. Macrophage expressed F13a1, Mrc1 and Apoe as M2 macrophage markers, while ACM4 expressed genes associated with pro-inflammatory microglia, such as cytokine release regulator Sp100 and Parp14. Interestingly, ACM3 exists at both GS and SE stages and the highly expressed gene in this cluster is enriched in nuclear division, DNA replication and chromosome segregation. Proliferation markers such as Pcna and Pola1 are the marker genes of ACM3. Therefore, ACM3 revealed an abnormal proliferative microglia type (Figure 6 and Table 3).

# Discussion 
FCD can induce drug-resistant epilepsy. While patients with mTOR signaling upregulation is associated with the seizures, how patients with such mutation can cause the disease is still unknown. This project is trying to find out potential molecular targets for future therapy.

Transcriptional factors as major gene expression regulators during neurodevelopment, may participate in mTOR induced targets critical for epilepsy development. Here we listed a series of TFs specifically expressed at SE stage. The expression of several targets specifically expressed in the neuron clusters were further verified pseudo-analysis. By reviewing the paper we identified several TFs related to neurodevelopment. For example, Id2 and Rai14 is identified to regulate axonal organization and dendritic spine dynamics, while Klf5 was reported to play roles in neurotransmitter signaling production and degeneration (Huang Z et al., 2019; Kim SJ et al., 2022, Wang Y st al., 2022). Thus, it was inferred that the re-organization of the synapse development by these TFs may be associated with the transition from GS to SE stage. 

Further sn-RNA data re-clustering identified an abnormal proliferative microglia cluster at the epilepsy stage. Thus, microglia exhibited a process of self-renewal in disease states. Previous studies revealed the activity dependent gene causes microglia proliferation after spinal cord injury (Tan et al., 2022). It was inferred that the activity stressed-induced signaling from ExN interacts with microglia to induce microgliosis after disease. 

This study has several limitations to be improved. Firstly, our animal epilepsy model is constructed through Pten knockout induced mTOR activation. Such model may be less applicable for FCD type I that is less influenced by mTOR. Other epilepsy mice models should be considered to find out shared molecular mechanism among models. Secondly, verification of biological experments such as RT-qPCR can be included to test the correctness of our data analysis. 

# Figure and tables
All figures and tables will be presented from Qingyuan's github: ee282. Figures will be separated into different parts:
Figure 1 part: Results of quality control;

Figure 2 part: Creating a single cell atlas of mice neocortex at different stages of epilepsy

Figure 3 part: Reclustering of ExN cell type and Argur analysis

Figure 4 part:  differential gene expression analysis from ExN

Figure 5 part: Pseudo-time analysis of layer II/III ExN

Figure 6 part: A self-renewal in microglia was identified

Table 1/2: Highly expressed genes from SE or GS compared with control

Table 3: Highly expressed genes of each microglia subcluster by DEG analysis

# Reference
Reference
Chen, Q.L. et al. (2020) ‘Bioinformatic analysis identifies key transcriptome signatures in temporal lobe epilepsy’, CNS Neuroscience and Therapeutics, 26(12), pp. 1266–1277. doi:10.1111/cns.13470.

Marsan, E. and Baulac, S. (2018) ‘Review: Mechanistic target of rapamycin (mTOR) pathway, focal cortical dysplasia and epilepsy’, Neuropathology and Applied Neurobiology, 44(1), pp. 6–17. doi:10.1111/nan.12463..

Pfisterer, U. et al. (2020) ‘Identification of epilepsy-associated neuronal subtypes and gene expression underlying epileptogenesis’, Nature Communications, 11(1). doi:10.1038/s41467-020-18752-7.

Singh, V.P. (2018) ‘Focal cortical dysplasia’, Neurology India, 66(6), pp. 1601–1602. doi:10.4103/0028-3886.246237.

Sisodiya, S.M. et al. (2009) ‘Focal cortical dysplasia type II: biological features and clinical perspectives’, The Lancet Neurology, 8(9), pp. 830–843. doi:10.1016/S1474-4422(09)70201-7.

Tan, W. et al. (2022) ‘Distinct phases of adult microglia proliferation: a Myc-mediated early phase and a Tnfaip3-mediated late phase’, Cell Discovery, 8(1). doi:10.1038/s41421-022-00377-3.

Wang Y, Cui Y, Liu J, Song Q, Cao M, Hou Y, Zhang X, Wang P. Krüppel-like factor 5 accelerates the pathogenesis of Alzheimer's disease via BACE1-mediated APP processing. Alzheimers Res Ther. 2022 Jul 26;14(1):103. doi: 10.1186/s13195-022-01050-3IF: 7.9 Q1 . PMID: 35883144; PMCID: PMC9316766.

Kim SJ, Woo Y, Kim HJ, Goo BS, Nhung TTM, Lee SA, Suh BK, Mun DJ, Kim JH, Park SK. Retinoic acid-induced protein 14 controls dendritic spine dynamics associated with depressive-like behaviors. Elife. 2022 Apr 25;11:e77755. doi: 10.7554/eLife.77755IF: 6.4 Q1 . PMID: 35467532; PMCID: PMC9068211.

Huang Z, Liu J, Jin J, et al. Inhibitor of DNA binding 2 promotes axonal growth through upregulation of Neurogenin2. Experimental Neurology. 2019 Oct;320:112966. DOI: 10.1016/j.expneurol.2019.112966. PMID: 31145898.