#######################################################################################################
# 2022-07
# RNA-Seq Project
# By Xiong Ewan
#######################################################################################################

# fastqc 
fastqc -t 24 -o /opt/test/xiongyifan/rnaseq/qc -q /opt/test/xiongyifan/rnaseq/rawdata/test_*.gz &

# cutadapt
cutadapt -j 6 --times 1  -e 0.1  -O 3  --quality-cutoff 25  -m 55 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-o /opt/test/xiongyifan/rnaseq/cleandata/test_R1_cutadapt.fq.gz \
-p /opt/test/xiongyifan/rnaseq/cleandata/test_R2_cutadapt.fq.gz \
/opt/test/xiongyifan/rnaseq/rawdata/test_R1.fq.gz \
/opt/test/xiongyifan/rnaseq/rawdata/test_R2.fq.gz > /opt/test/xiongyifan/rnaseq/cleandata/test_cutadapt.log 2>&1 &


# fastp

cat id |while read id;
do
fastp --thread 64 --detect_adapter_for_pe \
--overrepresentation_analysis \
--trim_front1 2 --trim_front2 2 \
--cut_front --cut_tail --cut_window_size 3 \
--cut_mean_quality 30 \
-i ${id}_1.fastq.gz -I ${id}_2.fastq.gz \
-o /opt/test/xiongyifan/rnaseq/public_data_test/trim/${id}_1_trim.fastq.gz \
-O /opt/test/xiongyifan/rnaseq/public_data_test/trim/${id}_2_trim.fastq.gz \
-j /opt/test/xiongyifan/rnaseq/public_data_test/trim/trim_report/${id} \
-h /opt/test/xiongyifan/rnaseq/public_data_test/trim/trim_report/${id} \
-R "Fastp Trimming Report"
done





###############################################################################################
# build index
###############################################################################################

# bowtie index
cd /home/menghaowei/ngs_course/reference/bowtie_index
bowtie-build -t 6 ref_hg38.fa  ref_hg38.fa > bt_build.log 2>&1 & 

# bowtie2 index
cd /home/menghaowei/ngs_course/reference/bowtie2_index
bowtie2-build -t 6 ref_hg38.fa ref_hg38.fa > bt2_build.log 2>&1 & 

# hisat2
hisat2-build ref_hg38.fa ref_hg38.fa > hisat2_build.log 2>&1 &

# STAR
STAR --runThreadN 96 --runMode genomeGenerate \
--genomeDir /opt/test/xiongyifan/reference/star_ref \
--genomeFastaFiles /opt/test/xiongyifan/reference/seq_hg38.fa \
--sjdbGTFfile /opt/test/xiongyifan/reference/hg38_ensembl.gtf \
--sjdbOverhang 150 & 

#   /opt/test/xiongyifan/reference/annotation/hg38_ensemble.gtf(no_chr)
#   /opt/test/xiongyifan/reference/hg38_chr.gtf
#    




###############################################################################################
# alignment
###############################################################################################
tophat2 -o ./test_tophat2 -p 6 \
-G /home/menghaowei/ngs_course/reference/gtf/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf \
/home/menghaowei/ngs_course/reference/bowtie2_index/ref_hg38.fa \
./fix.fastq/test_R1_cutadapt.fq.gz \
./fix.fastq/test_R2_cutadapt.fq.gz > test_tophat2/test_tophat2.log 2>&1 & 


hisat2 -p 12 \
-x /opt/test/xiongyifan/reference/hisat2_ref/ref_hg38 \
-1 /opt/test/xiongyifan/rnaseq/cleandata/test_R1_cutadapt.fq.gz \
-2 /opt/test/xiongyifan/rnaseq/cleandata/test_R2_cutadapt.fq.gz \
-S ./bam/test_hisat2.sam > ./bam/test_hisat2.log 2>&1 &


cat id |while read id;
do \
STAR \
--genomeDir /opt/test/xiongyifan/reference/star_ref \
--runThreadN 128 \
--readFilesIn ${id}_1_trim.fastq.gz ${id}_2_trim.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /opt/test/xiongyifan/rnaseq/public_data_test/align/bam/${id}_ \
--outSAMtype BAM Unsorted \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--outFilterIntronMotifs RemoveNoncanonical > /opt/test/xiongyifan/rnaseq/public_data_test/align/${id}_STAR.log 2>&1 
done

# STAR alignment & STAR-Fusion input Chimeric.out.junction
#   STAR 2-pass mode

STAR --genomeDir ${star_index_dir} \
          --runThreadN 128 \                                                                          
          --readFilesIn ${left_fq_filename} ${right_fq_filename} \  
          --readFilesCommand zcat \
          --outFileNamePrefix ${id}\                                                                    
          --outReadsUnmapped None \
          --twopassMode Basic \ ###
          --outSAMtype BAM Unsorted \
          --outSAMstrandField intronMotif \  # include for potential use with StringTie for assembly
          --outSAMattributes All \
          --outFilterIntronMotifs RemoveNoncanonical
          --outSAMunmapped Within \
          ############################
          ### STAR-Fusion relevant ### 
          ############################
          --chimSegmentMin 12 \  # ** essential to invoke chimeric read detection & reporting **
          --chimJunctionOverhangMin 8 \
          --chimOutJunctionFormat 1 \   # **essential** includes required metadata in Chimeric.junction.out file.
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
          --alignIntronMax 100000 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
          --outSAMattrRGline ID:GRPundef \
          --chimMultimapScoreRange 3 \
          --chimScoreJunctionNonGTAG -4 \
          --chimMultimapNmax 20 \
          --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 \
          --peOverlapMMp 0.1 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30


### BAM QC by qualimap
cat id |while read id;
do
qualimap rnaseq -bam /opt/test/xiongyifan/rnaseq/public_data_test/align/${id}_out.bam \
-gtf /opt/test/xiongyifan/reference/hg38_ensembl.gtf \
-pe \
-outdir /opt/test/xiongyifan/rnaseq/public_data_test/qc/bamqc/stroke_1 2>&1
done

##samtools.SAMFormatException: Error parsing SAM header. @RG line missing SM tag. Line:@RG

samtools view -H your.bam | grep -v "^@RG" | samtools reheader - your.bam > your.new.bam
samtools view -h ctr_1_test.bam|grep -v "^@RG" |samtools reheader - ctr_1_test.bam >t1.bam

samtools view -h ctr_1_test.bam|grep -v "^@RG" > t1.bam

RG:Z:GRPundef

# download some resource
## SNP
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

## GTF


# make exon 
hisat2_extract_exons.py hg38_refseq.gtf > hg38_refseq.exon &

# make splice site
hisat2_extract_splice_sites.py hg38_refseq.gtf > hg38_refseq.ss &

# make snp and haplotype
hisat2_extract_snps_haplotypes_UCSC.py ref_hg38.fa snp151Common.txt snp151Common &

# build index
hisat2-build -p 6 --snp snp151Common.snp --haplotype snp151Common.haplotype --exon hg38_refseq.exon  --ss hg38_refseq.ss ref_hg38.fa ref_hg38.fa.snp_gtf > hisat2_build.log 2>&1 & 

# mapping
hisat2 -p 6 \
-x /opt/test/xiongyifan/reference/hisat2_ref/ref_hg38 \
-1 /opt/test/xiongyifan/rnaseq/cleandata/test_R1_cutadapt.fq.gz \
-2 /opt/test/xiongyifan/rnaseq/cleandata/test_R2_cutadapt.fq.gz \
-S /opt/test/xiongyifan/rnaseq/alignment/test_hisat2.sam > /opt/test/xiongyifan/rnaseq/alignment/test_hisat2.log 2>&1 &


# sam 2 sorted bam

samtools sort -@ 64 -O bam -o output_name.bam 

###############################################################################################
# count by HTSeq and featureCount
###############################################################################################
# install 
conda install htseq
conda install subread

mamba install htseq
mamba install subread


# feature count
featureCounts -t exon -g gene_id \
-Q 10 --primary -s 0 -p -T 64 \
-a /opt/test/xiongyifan/reference/hg38_ensembl.gtf \
-o /opt/test/xiongyifan/rnaseq/counts/test_count.featureCounts \
/opt/test/xiongyifan/rnaseq/alignment/test_sorted.bam > /opt/test/xiongyifan/rnaseq/counts/test_count.featureCounts.log  2>&1 & 


 

# counts qc
qualimap counts -d config.txt -c -s human -outdir ./countqc_out &


-c,--compare             Perform comparison of conditions. Currently 2 maximum
                         is possible.
-d,--data <arg>          File describing the input data. Format of the file is
                         a 4 column tab-delimited table.
                         Column 1: sample name
                         Column 2: condition of the sample
                         Column 3: path to the counts data for the sample
                         Column 4: index of the column with counts




#novel transcripts 

stringtie \
-G /opt/test/xiongyifan/reference/hg38_ensembl.gtf \
-o /opt/test/xiongyifan/rnaseq/public_data_test/novel_transcript/${id}.gtf \
-p 128 \
/opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/name_sort.bam

stringtie --merge -p 96 -o all_sample.gtf *gtf &


# gffcompare

gffcompare -R -r /opt/test/xiongyifan/reference/hg38_ensembl.gtf -o outprefix stringtie_merge.gff



https://mirrors.tuna.tsinghua.edu.cn/bioconductor/3.13/





# Alternative splicing rMATS
conda activate rmats
rmats.py --b1 ctr.txt --b2 stroke.txt \
--gtf /opt/test/xiongyifan/reference/hg38_ensembl.gtf \
-t paired --readLength 100 --variable-read-length \
--cstat 0.0001 \
--nthread 64 \
--od ./out_put --tmp ./ &




~/anaconda3/envs/rmats/lib

/root/anaconda3/envs/rmats/rMATS/rMATS_C/rMATSexe: error while loading shared libraries: libcblas.so.3: cannot open shared object file: No such file or directory
从eca的lib cp libcblas.so.3 到rmats环境的lib
从eca的lib cp libcblas.so.3 到rmats环境的lib

rmats.py --b1 ctr.txt --b2 stroke.txt --gtf /opt/test/xiongyifan/reference/hg38_ensembl.gtf -t paired --readLength 147 --variable-read-length --cstat 0.0001 --nthread 64 --od ./out_put --tmp ./tmp_output &




conda activate plotrmats


convert a GTF to a GFF3 with gffread: gffread --keep-genes ./some_file.gtf -o ./some_file.gff3

# test
rmats2sashimiplot \
--b1 ./sample_1_replicate_1.bam,./sample_1_replicate_2.bam,./sample_1_replicate_3.bam \
--b2 ./sample_2_replicate_1.bam,./sample_2_replicate_2.bam,./sample_2_replicate_3.bam \
--l1 SampleOne --l2 SampleTwo --exon_s 1 --intron_s 5 -o plot_out -t SE -e ./SE.MATS.JC.txt &

# 
rmats2sashimiplot \
--b1 /opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/ctr_1_sort.bam,/opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/ctr_2_sort.bam,/opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/ctr_3_sort.bam \
--b2 /opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/stroke_1_sort.bam,/opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/stroke_2_sort.bam,/opt/test/xiongyifan/rnaseq/public_data_test/align/bam/sort_bam/stroke_3_sort.bam \
--l1 ctr --l2 stroke -o ./ctr_stroke_plot \
-t SE --exon_s 1 --intron_s 5 \
-e ./1.txt &






/opt/test/xiongyifan/reference/annotation/hg38_ensm.gtf

# hg38基因组_chr
/opt/test/genome_annotation/hg38_ensemble.gtf



#统计测序深度、覆盖度
#hg38 序列染色体大小
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
/opt/test/genome_annotation/hg38_chrom_sizes.txt


grep -vi 'random' /opt/test/genome_annotation/hg38_chrom_sizes.txt  | grep -vi un | while read i;do 
    conf=($i)
    chr=${conf[0]}
    len=${conf[1]}
    samtools depth -r $chr test.sorted.bam | \
    awk -v chr=$chr -v len=$len \
    '{sum+=$3} END { print chr,"-Average = ",sum/len}';done

# chr1 -Average =  1.01466
# chr2 -Average =  1.06905


awk '{sum+=$3} END { print chr,"-Average = ",sum/len}'



awk '{sum += $3};END {print sum}' chr1_depth.txt

for i in `awk '{$1}' hg38_chrom_sizes.txt`;while read i;do \
awk '{if($1==$i);sum+=$3};END {print sum}' chr1_depth.txt;done





/opt/test/bamdst/bamdst/bamdst -p \
/opt/test/pacbio/output/3_samtools_bam_sort/chr6.bed \
-o ./ /opt/test/pacbio/output/3_samtools_bam_sort/test.sorted.bam &

# bamdst 
# 配置bed 文件
# 能统计给定区域的深度、覆盖度。全染色体范围不适用


#融合基因
#STAR-fusion
# Visualization by R-package chimeraviz

# build STAR-fusion index

/root/anaconda3/envs/starfusion/bin/prep_genome_lib.pl \
                         --genome_fa ref_genome.fa \
                         --gtf gencode.*.annotation.gtf \
                         --fusion_annot_lib fusion_lib.*.dat.gz \
                         --annot_filter_rule AnnotFilterRule.pm \ #需额外下载3.5k
                         --pfam_db current \
                         --dfam_db human \
                         --human_gencode_filter #


/opt/test/xiongyifan/rnaseq/public_data_test/fusion/GRCh38_gencode_v37_CTAT_
STX16-NPEPL1--RAE1    4                  24                 INCL_NON_REF_SPLICE  STX16-NPEPL1^ENSG00000254995.4  chr20:57227143:+  RAE1^ENSG00000101146.8           chr20:55929088:+  YES_LDAS            5571.0306   GT              1.9899            AG               1.9656             INTRACHROMOSOMAL[chr20:1.27Mb]
lib_Mar012021.source

~/anaconda3/envs/starfusion/lib/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
--genome_fa GRCh38.primary_assembly.genome.fa \
--gtf gencode.v37.annotation.gtf \
--fusion_annot_lib fusion_lib.Mar2021.dat.gz\
--annot_filter_rule AnnotFilterRule.pm \
--pfam_db current \
--dfam_db human &



STAR-Fusion --genome_lib_dir /opt/test/xiongyifan/rnaseq/public_data_test/fusion/CTAT_resource_lib/ctat_genome_lib_build_dir \
             --left_fq /opt/test/xiongyifan/rnaseq/public_data_test/trim/ctr_1_1_trim.fastq \
             --right_fq /opt/test/xiongyifan/rnaseq/public_data_test/trim/ctr_1_2_trim.fastq \
             --output_dir fusion_out \
             --CPU 64



fq_site = /opt/test/xiongyifan/rnaseq/public_data_test/trim


STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib \
             -J Chimeric.out.junction \
             --output_dir star_fusion_outdir


${id}_1_trim.fastq


###########################################

# 染色体名不含chr，在染色体名前加上chr

###########################################

# gtf 文件
sed -i '/^#/d' GRCh38.gtf (去除#开头的注释信息)

sed -i '/^MT/d' GRCh38.gtf (去除MT染色体)

sed -i '/^KI/d' GRCh38.gtf （去除KI、GL开头的片段）

awk '{print "chr"$0}' GRCh38.gtf > GRCh38_chr.gtf  (加上chr)


# fasta 文件

sed 's/>/>chr/g' GRCh38.fa > GRCh38_chr.fa


# bam 文件
samtools view -H test.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - test.bam > test.CHR.bam


samtools view -H stroke_1_chr_spn.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - stroke_1_chr_spn.bam > stroke_1_change_chr.bam

###########################################

SNP/Indel calling by GATK

###########################################

conda activate cfDNApipe2.0

# STAR two pass mode mapping
STAR \
--genomeDir /opt/test/xiongyifan/reference/star_ref \
--runThreadN 128 \ 
--readFilesIn stroke_1_1_trim.fastq stroke_1_2_trim.fastq \
--outFileNamePrefix /opt/test/xiongyifan/rnaseq/public_data_test/snp_indel/stroke_1 \
--outSAMtype BAM Unsorted \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--twopassMode Basic \
--outFilterIntronMotifs RemoveNoncanonical > /opt/test/xiongyifan/rnaseq/public_data_test/snp_indel/stroke_1_STAR.log

# samtools sort bam 

samtools sort -@ 64 -O bam -o output_name.bam 


# gatk markduplicates remove bam duplicates
gatk --java-options "-Xmx6G" MarkDuplicates --REMOVE_DUPLICATES true \
--INPUT stroke_1_sort.bam  \
--METRICS_FILE stroke_1_rmdup.txt \
--OUTPUT stroke_1_rmdup.bam 2>&1 &


# create a .dict file form reference fastq
gatk CreateSequenceDictionary -R hg38.fa


# gatk SplitNCigarReads
gatk SplitNCigarReads \
-R /opt/test/xiongyifan/reference/seq_hg38.fa \
-I stroke_1_rmdup.bam \
-O stroke_1_splitn.bam &

hg38_ref = /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta


### Add @RG

gatk --java-options "-Xmx6G" AddOrReplaceReadGroups \
-I  stroke_1_chr_splitn.bam  \
-O  stroke_1_splitn_RG.bam \
-RGLB stroke_1 \
-RGPL ILLUMINA \
-RGPU stroke_1 \
-RGSM stroke_1


gatk IndexFeatureFile -I xxx.vcf.gz

#############
gatk BaseRecalibrator  \
-I stroke_1_splitn_RG.bam\
-R /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta  \
--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/1000G_omni2.5.hg38.vcf \
--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/1000G_phase1.snps.high_confidence.hg38.vcf  \
--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/hapmap_3.3.hg38.vcf \
--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/dbsnp_146.hg38.vcf \
--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-O stroke_1_recal.table

gatk ApplyBQSR \
--bqsr-recal-file stroke_1_recal.table \
-R /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta \
-I stroke_1_splitn_RG.bam \
-O stroke_1_recal.bam


# Variant Calling by HaplotypeCaller
gatk --java-options "-Xmx6g" HaplotypeCaller \
-I stroke_1_recal.bam \
-O stroke_1.vcf \
-R /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta \
--dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20
#-stand_call_conf相当于一个可信度打分，转录的默认是20，全基因组会考虑用30


# Variant filtering
gatk VariantFiltration \
-R /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta \
-V stroke_1.vcf \
-O stroke_1_filted.vcf \
-window 35 \
-cluster 3 \
-filter "FS > 30.0 || QD < 2.0" \
--filter-name FSQD 2>stroke_1_variantFilter.log



# variants annotation
/opt/test/pacbio/annovar/table_annovar.pl \
stroke_1_filted.vcf \
/opt/test/pacbio/annovar/humandb \
--out stroke_1_annovar.vcf \
-buildver hg38 -otherinfo -nastring . \
-protocol refGeneName,refGene,mim2gene,Gencode,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp42a \
--operation r,g,r,r,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f \
--vcfinput  --thread 30 --remove --intronhgvs 300 --argument '--colsWanted 5,--transcript_function,,--colsWanted 5,,,,,,,,,,,,,,,,' &
