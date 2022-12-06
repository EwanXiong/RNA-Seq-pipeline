rm(list = ls())
options(warn = -1)
suppressMessages(suppressWarnings(library(GetoptLong)))
GetoptLong(
    "outputdir=s", "output file site",
    "input=s","input file,must absolute path"
)
workdir = outputdir
setwd(workdir)
input=input
count <- read.table(input,header = T,sep = "\t")
trim_count <- count[,c(1,6:dim(count)[2])]
write.table(trim_count,paste0(outputdir,"/all_sample_count.txt"),quote = F,row.names = F,sep = '\t')
