去掉featurecount结果中的下游分析无关列
Rscript /opt/test/xiongyifan/script_R/trim_count_matrix.r \
--outputdir /opt/tsinghua/NuoHe/mRNA_22.11.24/output/6_counts \
--input all_counts

##将转录组下游差异分析分装
####CONVERT用法
Rscript /opt/tsinghua/test_file/transform/Final/convert_ID.r \
--outputdir ./ \
--input /opt/tsinghua/NuoHe/mRNA_22.11.17/output/DEGs/all_counts.txt \
--type mm \
####DEG用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/DEG.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--DEG 
####GO用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/GO.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--GO 
####KEGG用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/KEGG.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--KEGG 
####GSEA用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/GSEA.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--GESA 
####ALL总体流程用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/Auto_Enrichment_Final.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--ALL 
####修改config文件即可，注意config文件中的":"前后需要空格











