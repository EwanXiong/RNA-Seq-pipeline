##########################################################
##               DEG_v2                  ##
##########################################################
# 修改时间：2022.11.22
# 增补内容：1.TPM value
#        2.sample correlation heatmap plot
#        3.sample PCA plot
# 修改内容: 1.volcano plot 修改绘图数据，补充 not-siginficant的点


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
    "DEG!", "DEG KEGG GO enrichment analysis, Default: FALSE",
    "config=s","Config file, Used to set various parameters and variables"
)
config_file=config
config = yaml.load_file(config_file)
input_1=config$input_1
input_2=config$input_2
output_dir=config$output_dir


workdir = output_dir
setwd(workdir)


count <- read.table(file=input_1,header = T,row.names=1,sep = "\t")
sample_df <-read.table(file=input_2,header = T,row.names = 1)
count<-count[,c("Length",rownames(sample_df))]
########### DEG #############
deg<-function(){
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
    table(res$padj < 0.05)        # filter
    res <- res[order(res$padj),]  #按照padj进行升序排列
    res_plot <- res # Volcano plot use
    resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)
    diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1) 
    
    write.table(diff_gene_deseq2,file= paste0(output_dir,'/DEG',"/DEGs.txt"), quote = F,sep = "\t")
    cat("DEG results save finish\n")
    cat("FPKM analyze start\n")
    count_filter<-(count)[apply(count,1,function(x){ all(x > 1) }),]
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
    ggsave(paste0(output_dir,'/Expression','/sample_FPKM_heatmap.pdf'),ph,width = 6,height = 8)
    cat("plot correlation heatmap\n")
    cor_fpkm <- all_fpkm
    cor <- as.data.frame(cor(all_fpkm))
    cor_plot <- pheatmap(cor,treeheight_col = 15, treeheight_row=25,width = 30,height = 100,silent = T)
    ggsave(paste0(output_dir,'/Expression','/sample_correlation_heatmap.pdf'),cor_plot,width = 8,height = 7)
    write.table(cor,file= paste0(output_dir,'/Expression','/sample_cor.txt'), quote = F,sep = "\t")
    cat("Correlation analyze end\n")
    pca <- fast.prcomp(cor_fpkm)
    pca.data <- data.frame(sample = rownames(sample_df),Treatment = sample_df$condition,pca$rotation)
    pca_plot <- ggscatter(pca.data,x = "PC1",y = "PC2", color = 'Treatment',ellipse = T)+ theme_base()
    ggsave(paste0(output_dir,'/Expression','/sample_PCA.pdf'),pca_plot,width = 8,height = 7)
    cat("PCA analyze end\n")
    cat("All analyze end\n")
    return(res_plot)
}
summary_plot<-function(res){
        cat("Volcano plot start\n")
        #res <-read.table(paste0(output_dir,'/DEG',"/DEGs.txt"),header = T,row.names = 1)
        dir.create(path = paste0(output_dir,"/Diff"))
        up_DEG <- subset(res, padj < 0.05 & log2FoldChange > 1)
        down_DEG <- subset(res, padj < 0.05 & log2FoldChange < -1)
        write.table(up_DEG, paste0(output_dir,'/Diff','/Up_genes.txt'), quote = F,sep = "\t")      
        write.table(down_DEG, paste0(output_dir,'/Diff',"/Down_genes.txt"), quote = F,sep = "\t")
        voldata <- res
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
        ggsave(paste0(output_dir,"/Diff","/volcano_plot.pdf"),pp, width = 350,height = 350,units = 'mm')
    cat("Volcano plot finish\n")
}

if (DEG){
    res <- deg()
    summary_plot(res)
}
