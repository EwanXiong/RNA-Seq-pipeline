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
    "GO!", "DEG KEGG GO enrichment analysis, Default: FALSE",
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
if (!is.null(universe)){
  background = read.table(universe, header=FALSE, sep="\t")
  universe = as.character(background[,1])
}
go <- function(){
    cat("Go analyze start,start select specis and loading DEG results\n")
    dir.create(path = paste0(output_dir,"/GO"))
    if (species=="Hs"){
        lapply(c('BP', 'CC', 'MF'), function(ont){
          if (type != "ENTREZID"){
              idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Hs.eg.db")
              IDs <- idss$ENTREZID
              type <- "ENTREZID"
            }
          ego <- enrichGO(IDs, 'org.Hs.eg.db', keyType = type, ont = ont, pvalueCutoff = 0.05, 
          pAdjustMethod = pAdjustMethod, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500, universe = universe, readable = TRUE, pool = FALSE)
          if (!is.null(ego)){
            try(write.csv(as.data.frame(ego),paste0(output_dir,"/GO","/GO_enrichment_",ont,".csv"),row.names = F)) 
            try(dotplot(ego, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/GO","/dotplot_",ont,".pdf"),.,width = 8,height = 10))
            try(dag<-plotGOgraph(ego))
            try(ggsave(paste0(output_dir,"/GO","/GOgraph_",ont,".pdf"),dag,width = 8,height = 10))
            try(barplot(ego, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/GO","/barplot_",ont,".pdf"),.,width = 8,height = 10))
            ego_simp <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
            write.csv(as.data.frame(ego_simp),paste0(output_dir,"/GO","/GO_enrichment_simplified_",ont,".csv"),row.names = F)
            return(paste0("Genes identified for GO enrichment of sub-ontology ",ont,"."))
          }
            else{
            return(paste0("No genes are enriched for sub-ontology ",ont,"."))
          }
        }
        )
    }
    else if (species=="Mm"){
        lapply(c('BP', 'CC', 'MF'), function(ont){
          if (type != "ENTREZID"){
              idss <- bitr(IDs, fromType=type, toType="ENTREZID", OrgDb="org.Mm.eg.db")
              IDs <- idss$ENTREZID
              type <- "ENTREZID"
            }
          ego <- enrichGO(IDs, 'org.Mm.eg.db', keyType = type, ont = ont, pvalueCutoff = 0.05, 
          pAdjustMethod = pAdjustMethod, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500, universe = universe, readable = TRUE, pool = FALSE)
          if (!is.null(ego)){
            try(write.csv(as.data.frame(ego),paste0(output_dir,"/GO","/GO_enrichment_",ont,".csv"),row.names = F)) 
            try(dotplot(ego, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/GO","/dotplot_",ont,".pdf"),.,width = 8,height = 10))
            try(dag<-plotGOgraph(ego))
            try(ggsave(paste0(output_dir,"/GO","/GOgraph_",ont,".pdf"),dag,width = 8,height = 10))
            try(barplot(ego, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/GO","/barplot_",ont,".pdf"),.,width = 8,height = 10))
            ego_simp <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
            write.csv(as.data.frame(ego_simp),paste0(output_dir,"/GO","/GO_enrichment_simplified_",ont,".csv"),row.names = F)
            return(paste0("Genes identified for GO enrichment of sub-ontology ",ont,"."))
          }
            else{
            return(paste0("No genes are enriched for sub-ontology ",ont,"."))
          }
        }
        )
    }
    cat("GO analyze finish\n")
}

if (GO){
    go()
}
