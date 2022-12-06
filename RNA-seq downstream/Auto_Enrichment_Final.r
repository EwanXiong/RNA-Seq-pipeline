##########################################################
#           Auto_DEGs_Enrichment              #
##########################################################
# Analysis including:
# 1. DESeq2 
# 2. FPKM & TPM value
# 3. DEGs FPKM heatmap
# 4. Sample correlation heatmap
# 5. Sample PCA analyze
# 6. Volcano plot
# 7. GO , KEGG & GSEA enrichmemnt
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
suppressPackageStartupMessages(library(gmodels))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(org.Mm.eg.db)))
suppressMessages(suppressWarnings(library(enrichplot)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(yaml)))
GetoptLong(
    "ALL!", "DEG KEGG GO enrichment analysis, Default: FALSE",
    "config=s","Config file, Used to set various parameters and variables"
)

config_file=config
config = yaml.load_file(config_file)
input_1=config$input_1
input_2=config$input_2
output_dir=config$output_dir
species=config$species
type =config$type
pAdjustMethod = config$pAdjustMethod
universe=config$universe
showCategory = config$showCategory
KOannoation=config$KOannoation
GOannoation=config$GOannoation
TB=config$TB

workdir = output_dir
setwd(workdir)

count <- read.table(file=input_1,header = T,row.names=1,sep = "\t")
sample_df <-read.table(file=input_2,header = T,row.names = 1)
count<-count[,c("Length",rownames(sample_df))]
########### DEG #############
DEG<-function(){
    cat("DEG analyze start\n")
    dir.create(path = paste0(output_dir,"/DEG"))
    dir.create(path = paste0(output_dir,"/Expression"))
    dir.create(path = paste0(output_dir,"/Diff"))
    deg_count<-(count[,-1])[apply(count[,-1],1,function(x){ all(x > 1) }),]
    if(all(rownames(sample_df) == colnames(deg_count))){cat ("Sample ID is TRUE \n")}
    dds <- DESeqDataSetFromMatrix(countData = deg_count, colData = sample_df, design = ~condition)
    ### sample cluster
    vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
    sampleDists <- dist(t(assay(vsd))) 
    res1 <- hcut(sampleDists, k = length(unique(sample_df$condition)), stand = FALSE,hc_method ="average" ) 
    p<-fviz_dend(res1,rect_fill = T,cex = 1,color_labels_by_k=T,horiz=F)
    ggsave(paste0(output_dir,'/DEG','/sample_cluster.pdf'),p,width = 250,height = 200,units = 'mm')
    cat("Sample cluster plot finish\n")
    # normalization 
    dds <- estimateSizeFactors(dds)
    # dispersion
    dds <- estimateDispersions(dds)
    # plot dispersion
    pdf(paste0(output_dir,'/DEG','/sample_dispersion.pdf'), width = 14, height = 8)
    plotDispEsts(dds, ymin = 1e-4)
    dev.off()
    cat("sample dispersion plot finish\n")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    dds1 <- DESeq(dds)    # 数据标准化
    resultsNames(dds1)    # 查看结果的名称
    dds1$condition        #默认后者的处理组比前面的对照组
    res <- results(dds1)  
    summary(res)          #看一下结果的概要信息，adjusted p值默认 < 0.1
    table(res$pvalue < 0.05)        # filter
    res <- res[order(res$pvalue),]  #按照padj进行升序排列
    res_plot <- res # Volcano plot use
    resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)
    diff_gene_deseq2 <-subset(res, pvalue < 0.05 & abs(log2FoldChange) > 1)
    write.table(res,file= paste0(output_dir,'/DEG',"/DEGs_no_filter.txt"), quote = F,sep = "\t")# 不进行阈值筛选的DESeq2结果
    write.table(diff_gene_deseq2,file= paste0(output_dir,'/DEG',"/DEGs.txt"), quote = F,sep = "\t")
    cat("DEG results save finish\n")
    cat("FPKM analyze start\n")
    count_filter<-(count)[apply(count,1,function(x){ all(x >= 1) }),]
    count_df<-count_filter[,-1]
    countToFpkm <- function(counts, effLen)
    {
      N <- colSums(counts)
      exp( log(counts) + log(1e9) - log(effLen) - log(N) )
    }
    all_fpkm <- countToFpkm(count_df, count_filter$Length)
    write.table(all_fpkm,file= paste0(output_dir,'/Expression','/all_FPKM.txt'), quote = F,sep = "\t")
    cat("FPKM value save\n")
    countToTpm <- function(counts, effLen)
    {   
      rate <- log(counts) - log(effLen)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }   
    all_tpm <- countToTpm(count_df, count_filter$Length)
    write.table(all_tpm,file= paste0(output_dir,'/Expression','/all_TPM.txt'), quote = F,sep = "\t")
    cat("TPM value save\n")
    cat("FPKM pheatmap plot\n")
    deg_fpkm <-all_fpkm[rownames(all_fpkm)%in%rownames(diff_gene_deseq2),]
    deg_fpkm_log2 <- log2(deg_fpkm)
    ph <- pheatmap(deg_fpkm_log2+1,show_rownames = F,color = colorRampPalette(c("skyblue","white","red"))(256),
         treeheight_col = 15, treeheight_row=25,width = 30,height = 100,silent = T)
    ggsave(paste0(output_dir,'/Expression','/sample_FPKM_heatmap.pdf'),ph,width = 6,height = 8 )
    cat("plot correlation heatmap\n")
    cor_fpkm <- all_fpkm
    cor <- as.data.frame(cor(all_fpkm))
    cor_plot <- pheatmap(cor,treeheight_col = 15, treeheight_row=25,width = 30,height = 100,silent = T)
    ggsave(paste0(output_dir,'/Expression','/sample_correlation_heatmap.pdf'),cor_plot,width = 8,height = 7)
    write.table(cor,file= paste0(output_dir,'/Expression','/sample_cor.txt'), quote = F,sep = "\t")
    cat("Correlation analyze end\n")
    cat("PCA analyze\n")
    pca <- fast.prcomp(cor_fpkm)
    pca.data <- data.frame(sample = rownames(sample_df),Treatment = sample_df$condition,pca$rotation)
    pca_plot <- ggscatter(pca.data,x = "PC1",y = "PC2", color = 'Treatment',ellipse = T)+ theme_base()
    ggsave(paste0(output_dir,'/Expression','/sample_PCA.pdf'),pca_plot,width = 8,height = 7)
    cat("PCA analyze end\n")
    cat("All analyze end\n")
    return(res_plot)
}
summary_plot<-function(){
        cat("Volcano plot start\n")
        #res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1) res改由DEG()return出直接调用
        dir.create(path = paste0(output_dir,"/Diff"))
        up_DEG <- subset(res, padj < 0.05 & log2FoldChange > 1)
        down_DEG <- subset(res, padj < 0.05 & log2FoldChange < -1)
        write.table(up_DEG, paste0(output_dir,'/Diff','/Up_genes.txt'), quote = F,sep = "\t")      
        write.table(down_DEG, paste0(output_dir,'/Diff',"/Down_genes.txt"), quote = F,sep = "\t")
        voldata <- res[!is.na(res$padj),]
        voldata$threshold = factor(ifelse(voldata$padj < 0.05 & abs(voldata$log2FoldChange) >= 1,ifelse(voldata$log2FoldChange >= 1, 'Up-regulated','Down-regulated'),'not-siginficant'),levels = c('Up-regulated','Down-regulated','not-siginficant'))
        voldata$Gene <- rownames(voldata)
        vol_plot <- as.data.frame(voldata)
        pp<-ggplot(vol_plot, aes(
          x = log2FoldChange,
          y = -log10(padj),
          color = threshold
        )) +
          geom_point(alpha = 0.5, size = 3, shape = 19, fill = 'white') +
          scale_color_manual(values = c("#DC143C","#00008B","#808080")) + ##2f5688","#C0C0C0","#CC0000"
          geom_text_repel(
            #添加差异标签
            data = vol_plot[vol_plot$padj < 0.05 & abs(vol_plot$log2FoldChange) >2,],
            aes(label = Gene),
            size = 3,
            segment.color = "black",
            show.legend = FALSE
          ) + #添加关注的点的基因id
          theme_bw() +
          theme(legend.title = element_blank()) +
          xlim(-7.5, 7.5) +
          ylab('-log10(P.adjusted Value)') +
          xlab('log2(FoldChange)') +
          geom_vline(
            xintercept = c(-1,1),
            lty = 3,
            col = "black",
            lwd = 0.5
          ) + #添加横线
          geom_hline(
            yintercept = -log10(0.05),
            lty = 3,
            col = "black",
            lwd =0.5
          )
        ggsave(paste0(output_dir,"/Diff","/volcano_plot.pdf"),pp, width = 300,height = 200,units = 'mm')
    cat("Volcano plot finish\n")
}
universe=NULL
if (!is.null(universe)){
  background = read.table(universe, header=FALSE, sep="\t")
  universe = as.character(background[,1])
}
go <- function(){
    cat("Go analyze start,start select specis and loading DEG results\n")
    dir.create(path = paste0(output_dir,"/GO"))
    if (species=="Hs"){
        res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
        IDs <- as.character(rownames(res))
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
            try(plotGOgraph(ego)%>%ggsave(paste0(output_dir,"/GO","/GOgraph_",ont,".pdf"),.,width = 8,height = 10))
            #try(pdf(paste0(output_dir,"/GO","/GO_DAGplot_",ont,".pdf"),width = 10,height = 8)) 
            #try(plotGOgraph(ego))  
            #try(dev.off())  
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
        res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
        IDs <- as.character(rownames(res))
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
            #try(pdf(paste0(output_dir,"/GO","/GO_DAGplot_",ont,".pdf"),width = 10,height = 8)) 
            #try(plotGOgraph(ego))  
            #try(dev.off())  
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
    else {
        res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
        dir.create(path = paste0(output_dir,"/GO"))
        gene_list <- as.character(rownames(res))
        KOannotation <- read.delim(KOannoation, stringsAsFactors=FALSE)
        GOannotation <- read.delim(GOannoation, stringsAsFactors=FALSE)
        GOinfo <- read.delim(TB, stringsAsFactors=FALSE)
        GOannotation = split(GOannotation, with(GOannotation, level))
        for (ont in c('BP', 'CC', 'MF')){
            ego<-enricher(gene_list,TERM2GENE=GOannotation[[ont]][c(2,1)],TERM2NAME=GOinfo[1:2])
            if (!is.null(ego)){
            try(write.table(ego,file= paste0(output_dir,"/GO","/GO_",ont,".txt"), quote = F,sep = "\t",row.names = F))
            try(dotplot(ego, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/GO","/dotplot_",ont,".pdf"),.,width = 8,height = 10))
            #try(plotGOgraph(ego)%>%ggsave(paste0(output_dir,"/GO","/GOgraph_",ont,".pdf"),.,width = 8,height = 10))
            try(pdf(paste0(output_dir,"/GO","/GO_DAGplot_",ont,".pdf"),width = 10,height = 8)) 
            try(plotGOgraph(ego))  
            try(dev.off()) 
            try(barplot(ego, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/GO","/barplot_",ont,".pdf"),.,width = 8,height = 10))
                }
        }
    }
    cat("GO analyze finish\n")
}
kegg <- function(){
    cat("KEGG analyze start,start select genelist from DEG results\n")
    dir.create(path = paste0(output_dir,"/KEGG"))
    res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
    IDs <- as.character(rownames(res))
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
        res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
        IDs <- as.character(rownames(res))
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
    else {
    res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
    IDs <- as.character(rownames(res))
    kegg_enrich = enricher(IDs,TERM2GENE=KOannoation[c(3,1)],TERM2NAME=KOannoation[c(3,4)])
    try(write.csv(kegg_enrich,paste0(output_dir,"/KEGG","/KEGG_enrichment.csv"),row.names =F))
    try(dotplot(kegg_enrich, showCategory=showCategory, x = "GeneRatio", color = "p.adjust", font.size = 12)%>%ggsave(paste0(output_dir,"/KEGG","/dotplot.pdf"),.,width = 8,height = 10))
    try(barplot(kegg_enrich, showCategory=showCategory)%>%ggsave(paste0(output_dir,"/KEGG","/barplot.pdf"),.,width = 8,height = 10))
        }
    cat("KEGG analyze finish\n")
}

gsea<-function(){
    cat("GESA analyze start,start select genelist from DEG results\n")
    dir.create(path = paste0(output_dir,"/GSEA"))
    res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)[2]
    IDs <- as.character(rownames(res))
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
        #try(gseaplot2(gse.KEGG,1:5,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_2.pdf"),.,width = 15,height = 10))
        #try(gseaplot2(gse.KEGG,1:10,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_3.pdf"),.,width = 15,height = 10))
        try(write.table(gse.KEGG@result,file= paste0(output_dir,"/GSEA","/GSEA_pathway_results.txt"), quote = F,sep = "\t"))
        try(write.table(gsea,file= paste0(output_dir,"/GSEA","/GSEA_enrichment_results.txt"), quote = F,sep = "\t"))
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
        #try(gseaplot2(gse.KEGG,1:5,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_2.pdf"),.,width = 15,height = 10))
        #try(gseaplot2(gse.KEGG,1:10,pvalue_table = T)%>%ggsave(paste0(output_dir,"/GSEA","/gseplot_3.pdf"),.,width = 15,height = 10))
        try(write.table(gse.KEGG@result,file= paste0(output_dir,"/GSEA","/GSEA_pathway_results.txt"), quote = F,sep = "\t"))
        try(write.table(gsea,file= paste0(output_dir,"/GSEA","/GSEA_enrichment_results.txt"), quote = F,sep = "\t"))
    }
    else cat("specise can not GSEA analyze!")
    cat("GSEA analyze finish\n")
}


if (ALL){
    DEG()
    res <- DEG()
    summary_plot()
    go()
    kegg()
    gsea()
}
