rm(list = ls())

options(warn = -1)
suppressMessages(suppressWarnings(library(GetoptLong)))
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(org.Mm.eg.db)))

GetoptLong(
    "outputdir=s", "output file site",
    "input=s","input file,must absolute path",
    "type=s","hsa or mm"
)
workdir = outputdir
setwd(workdir)
input=input
count <- read.table(input,header = T,sep = "\t")
count<-rename(count, "Geneid"= "ENSEMBL")
if (type == "hsa"){
    change_gene_ID = bitr(count$ENSEMBL, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
    change_gene_ID=change_gene_ID[!(duplicated(change_gene_ID$SYMBOL)),]
    new_count <- merge(change_gene_ID,count,"ENSEMBL")
    write.table(new_count,paste0(outputdir,"/ENS_SYM_count.txt"),quote = F,row.names = F,sep = '\t')
    write.table(new_count[-1],paste0(outputdir,"/SYM_count.txt"),quote = F,row.names = F,sep = '\t')
}else if (type == "mm"){
    change_gene_ID = bitr(count$ENSEMBL, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
    change_gene_ID=change_gene_ID[!(duplicated(change_gene_ID$SYMBOL)),]
    new_count <- merge(change_gene_ID,count,"ENSEMBL")
    write.table(new_count,paste0(outputdir,"/ENS_SYM_count.txt"),quote = F,row.names = F,sep = '\t')
    write.table(new_count[-1],paste0(outputdir,"/SYM_count.txt"),quote = F,row.names = F,sep = '\t')
}else cat("Specise Type erro\n")
