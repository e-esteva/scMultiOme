library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(glue)
library(dplyr)
library(tidyverse)
library(celldex)
library(SingleR)
library(SingleCellExperiment)


######### adapted from: 
## I. https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1
## II. https://stuartlab.org/signac/articles/pbmc_vignette.html


# cellranger output:
args=commandArgs(TRUE)

cellranger_out=args[1]
message(glue('Cellranger Outs: {cellranger_out}'))

out_dir=args[2]
qc.dir=glue('{out_dir}/QC/')
dir.create(qc.dir)

message(glue('QC dir: {qc.dir}'))

dr.dir=glue('{out_dir}/DimReduction/')
message(glue('DimReductions dir: {dr.dir}'))
dir.create(dr.dir)

anno.dir=glue('{out_dir}/SingleR-Annotations/')
message(glue('Annotations dir: {anno.dir}'))
dir.create(anno.dir)

species = args[3]
if(tolower(species) %in% c('mm10','mm9','musmusculus','mouse')){
	species = 'mouse'
	

}else{
	if(tolower(species) %in% c('hg38','hg39','hsapiens','human')){
        	species	= 'human'

	}
}
message(glue('species: {species}'))

#################
setwd(cellranger_out)

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")


# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

metadata.file='per_barcode_metrics.csv'
metadata <- read.csv(
  file = metadata.file,
  header = TRUE,
  row.names = 1
)

# Create Seurat object
obj <- CreateSeuratObject(counts = rna_counts,meta.data=metadata)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
if(species == 'mouse'){
	genome='mm10'
	annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

	seqlevelsStyle(annotations) <- 'UCSC'
	genome(annotations) <- genome
}else{
	genome='hg38'
	annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
	seqlevelsStyle(annotations) <- 'UCSC'
	genome(annotations) <- genome

}


frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = genome,
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
obj[["ATAC"]] <- chrom_assay

setwd(qc.dir)

# QC
qc.pre=VlnPlot(obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

ggsave('QC-prefilter.pdf',qc.pre)


DefaultAssay(obj)='ATAC'
# compute nucleosome signal score per cell
obj <- NucleosomeSignal(object = obj)

# compute TSS enrichment score per cell
obj <- TSSEnrichment(object = obj, fast = FALSE)


# add blacklist ratio and fraction of reads in peaks
# atac_peak_region_fragments == peak_region_fragments
# atac_fragments == passed_filters
obj$pct_reads_in_peaks <- obj$atac_peak_region_fragments / obj$atac_fragments * 100
# blacklist:https://github.com/Boyle-Lab/Blacklist

if(genome == 'mm10'){

	mm10.blacklist=read.table('blacklists/mm10.blacklist.bed',sep='\t',header = F)
	names(mm10.blacklist)=c('chr','start','end')
	mm10.blacklist=try(makeGRangesFromDataFrame(mm10.blacklist))
	#print(head(mm10.blacklist))
	if(!('mm10.blacklist-granges.rds' %in% list.files('blacklists'))){
		saveRDS(mm10.blacklist,'blacklists/mm10.blacklist-granges.rds')
	}
	if(class(mm10.blacklist) == 'try-error'){
		mm10.blacklist=readRDS('blacklists/mm10.blacklist-granges.rds')

	}
	obj$blacklist_region_fragments=CountsInRegion(object = obj,regions = mm10.blacklist,assay = 'ATAC')
	obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$atac_peak_region_fragments
}else{
	hg38.blacklist=read.table('blacklists/hg38-blacklist.v2.bed',sep='\t',header = F)
        names(hg38.blacklist)=c('chr','start','end')
        hg38.blacklist=try(makeGRangesFromDataFrame(hg38.blacklist))

        if(!('hg38.blacklist-granges.rds' %in% list.files('blacklists/'))){
                saveRDS(hg38.blacklist,'blacklists/hg38.blacklist-granges.rds')
        }

        if(class(hg38.blacklist) == 'try-error'){
                hg38.blacklist=readRDS('blacklists/hg38.blacklist-granges.rds')

        }
        obj$blacklist_region_fragments=CountsInRegion(object = obj,regions = hg38.blacklist,assay = 'ATAC')
        obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$atac_peak_region_fragments



}
obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')

tss.plot=TSSPlot(obj, group.by = 'high.tss') + NoLegend()
ggsave('TSS-plot.pdf',tss.plot)

obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
if(length(unique(obj$nucleosome_group)) > 1){
	nuclesomeplot1=try(FragmentHistogram(object = obj, group.by = 'nucleosome_group'))
	ggsave('nucleosome-histogram-byGroup.pdf',nuclesomeplot1)
}

atac.qc=VlnPlot(
  object = obj,
  features = c('pct_reads_in_peaks', 'atac_peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave('ATAC-QC-panel.pdf',atac.qc)

# for better visibility:
atac.qc_adj=VlnPlot(
  object = obj,
  features = c('pct_reads_in_peaks', 'atac_peak_region_fragments'),
  pt.size = 0.1,
  ncol = 2
)
ggsave('QC-panel-adj-0.pdf',atac.qc_adj)

atac.qc_adj=VlnPlot(
  object = obj,
  features = c(
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 3
)
ggsave('QC-panel-adj-1.pdf',atac.qc_adj)

obj <- subset(
  x = obj,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.1 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

qc.summary_=data.frame('cells_pre_QC'=dim(rna_counts)[2],'cells_post_QC'=dim(obj)[2])
qc.summary_$`percent lost`=abs(qc.summary_$cells_post_QC/qc.summary_$cells_pre_QC-1)*100
write.csv(qc.summary_,'qc-summary.csv')

# for better visibility:
atac.qc_adj=VlnPlot(
  object = obj,
  features = c('pct_reads_in_peaks', 'atac_peak_region_fragments'),
  pt.size = 0.1,
  ncol = 2
)
ggsave('postQC-panel-adj-0.pdf',atac.qc_adj)

atac.qc_adj=VlnPlot(
  object = obj,
  features = c(
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 3
) 
ggsave('postQC-panel-adj-1.pdf',atac.qc_adj)


rna.qc_adj=VlnPlot(
  object = obj,
  features = c(
               "nFeature_RNA","nCount_RNA","percent.mt"),
  pt.size = 0.1,
  ncol = 3
) 
ggsave('postQC-panel-RNA.pdf',rna.qc_adj)


gene.activities <- GeneActivity(obj)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
obj[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(
  object = obj,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj$nCount_GeneActivity)
)


saveRDS(obj,'postQC-preNorm.rds')

setwd(dr.dir)

# Dim Redux independently
# RNA analysis
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

umap=DimPlot(obj,reduction='umap.rna')
ggsave('RNA-UMAP.pdf',umap)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

umap=DimPlot(obj,reduction='umap.atac')
ggsave('ATAC-UMAP.pdf',umap)

# We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities.
# We use this graph for UMAP visualization and clustering



obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p=DimPlot(obj,label=T,reduction='wnn.umap')
ggsave('WNN-UMAP.pdf',p)


saveRDS(obj,'scRNA-ATAC.rds')

####### singleR cellDex annotation
setwd(anno.dir)

source('getCellDex-annos.R')

obj=obj %>% NormalizeData(obj)

getCellDex(obj,species,assay='SCT')

