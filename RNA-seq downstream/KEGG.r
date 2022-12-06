rm(list = ls())

options(warn = -1)
suppressMessages(suppressWarnings(library(GetoptLong)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(factoextra)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggthemes))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(org.Mm.eg.db)))
suppressMessages(suppressWarnings(library(enrichplot)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(yaml)))

GetoptLong(
    "KEGG!", "DEG KEGG GO enrichment analysis, Default: FALSE",
    "config=s","Config file, Used to set various parameters and variables"
)
config_file=config
config = yaml.load_file(config_file)
input=config$input
output_dir=config$output_dir
species=config$species
type =config$type
pAdjustMethod = config$pAdjustMethod
universe=config$universe
showCategory = config$showCategory


workdir = output_dir
setwd(workdir)

ID=read.table(input,header = T)
IDs=as.character(ID[,1])
kegg <- function(){
    cat("KEGG analyze start,start select genelist from DEG results\n")
    dir.create(path = paste0(output_dir,"/KEGG"))
    if (species=="Hs"){
        if (type != "ENTREZID"){
          idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Hs.eg.db")
          IDs <- idss$ENTREZID
        }
        if (!is.null(universe) && type != "ENTREZID"){
            idsss <- bitr(universe, fromType=type, toType="ENTREZID", OrgDb="org.Hs.eg.db")
          universe <- idsss$ENTREZID
        }
    try(kegg_enrich <- enrichKEGG(gene=IDs, organism = 'hsa', keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = pAdjustMethod, qvalueCutoff = 1, universe = universe, minGSSize = 5, maxGSSize = 500, use_internal_data = FALSE))
    try(kk <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID"))
    try(write.csv(as.data.frame(kk),paste0(output_dir,"/KEGG","/KEGG_enrichment.csv"),row.names =F))
    try(dotplot(kegg_enrich, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/KEGG","/dotplot.pdf"),.,width = 8,height = 10))
    try(barplot(kegg_enrich, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/KEGG","/barplot.pdf"),.,width = 8,height = 10))
        }
    else if (species=="Mm"){
        if (type != "ENTREZID"){
          idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Mm.eg.db")
          IDs <- idss$ENTREZID
        }
        if (!is.null(universe) && type != "ENTREZID"){
          idsss <- bitr(universe, fromType=type, toType="ENTREZID", OrgDb="org.Mm.eg.db")
          universe <- idsss$ENTREZID
        }
    try(kegg_enrich <- enrichKEGG(gene=IDs, organism = 'mmu', keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = pAdjustMethod, qvalueCutoff = 1, universe = universe, minGSSize = 5, maxGSSize = 500, use_internal_data = FALSE))
    try(kk <- setReadable(kegg_enrich, OrgDb =org.Mm.eg.db, keyType="ENTREZID"))
    try(write.csv(as.data.frame(kk),paste0(output_dir,"/KEGG","/KEGG_enrichment.csv"),row.names =F))
    try(dotplot(kegg_enrich, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/KEGG","/dotplot.pdf"),.,width = 8,height = 10))
    try(barplot(kegg_enrich, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/KEGG","/barplot.pdf"),.,width = 8,height = 10))
        }
    cat("KEGG analyze finish\n")
}
if (KEGG){
    kegg()
}
