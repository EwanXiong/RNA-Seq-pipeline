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
    "GSEA!", "GSEA enrichment analysis, Default: FALSE",
    "config=s","Config file, Used to set various parameters and variables"
)

config_file=config
config = yaml.load_file(config_file)
input=config$input_gsea
output_dir=config$output_dir
species=config$species
type =config$type


workdir = output_dir
setwd(workdir)

res <-read.table(input,header = T,row.names = 1)[2]
IDs <- as.character(rownames(res)) 
gsea<-function(){
    cat("GESA analyze start,start select genelist from DEG results\n")
    dir.create(path = paste0(output_dir,"/GSEA"))
    if (species=="Hs"){
        if (type != "ENTREZID"){
          idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Hs.eg.db")
        }
        res$ENSEMBL<-rownames(res)
        genelist<-arrange(merge(idss,res)[2:3],desc(log2FoldChange))
        geneList<-genelist$log2FoldChange
        names(geneList)<-genelist$ENTREZID
        kegg_gmt <- read.gmt("/opt/tsinghua/software/GSEA/hsa/msigdb.v7.5.1.entrez.gmt")
        gsea <- GSEA(geneList,TERM2GENE = kegg_gmt)
        gse.KEGG <- gseKEGG(geneList, organism = "hsa",pvalueCutoff = 1,pAdjustMethod = "BH")
        try(dotplot(gsea)%>%ggsave(paste0(output_dir,"/GSEA","/dotplot.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:3,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_1.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:5,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_2.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:10,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_3.pdf"),.,width = 15,height = 10))
        try(write.table(gse.KEGG@result,file= paste0(output_dir,"/GSEA","/GSEA_KEGG_results.txt"), quote = F,sep = "\t"))
        try(write.table(gsea,file= paste0(output_dir,"/GSEA","/GSEA_results.txt"), quote = F,sep = "\t"))
    }
    if (species=="Mm"){
        if (type != "ENTREZID"){
          idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Mm.eg.db")
        }
        res$ENSEMBL<-rownames(res)
        genelist<-arrange(merge(idss,res)[2:3],desc(log2FoldChange))
        geneList<-genelist$log2FoldChange
        names(geneList)<-genelist$ENTREZID
        kegg_gmt <- read.gmt("/opt/tsinghua/software/GSEA/mmu/msigdb.v2022.1.Mm.entrez.gmt")
        gsea <- GSEA(geneList,TERM2GENE = kegg_gmt)
        gse.KEGG <- gseKEGG(geneList, organism = "mmu",pvalueCutoff = 1,pAdjustMethod = "BH")
        try(dotplot(gsea)%>%ggsave(paste0(output_dir,"/GSEA","/dotplot.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:3,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_1.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:5,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_2.pdf"),.,width = 15,height = 10))
        try(gseaplot2(gse.KEGG,1:10,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_3.pdf"),.,width = 15,height = 10))
        try(write.table(gse.KEGG@result,file= paste0(output_dir,"/GSEA","/GSEA_KEGG_results.txt"), quote = F,sep = "\t"))
        try(write.table(gsea,file= paste0(output_dir,"/GSEA","/GSEA_results.txt"), quote = F,sep = "\t"))
    }
    cat("GSEA analyze finish\n")
}

if (GSEA){
    gsea()
}
