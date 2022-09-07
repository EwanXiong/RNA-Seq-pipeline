####   RNA-seq shell pipeline   #####
# 2022-08-22
# Ewan Xiong 
#####################################

### load file site & var prepare ###
input=/opt/test/xiongyifan/rnaseq/auto_rna/input
output=/opt/test/xiongyifan/rnaseq/auto_rna/output
sample_name=($(cd $input; ls *gz|sed 's/_..fastq.gz//g'|uniq))

if [ ! -d "$output/1_qc" ]; then
	echo '-----------------------------------------------------------'
	echo "start time: `date` | procress 1_qc"
	echo '-----------------------------------------------------------'
	mkdir -p $output/1_qc
	sample=($(cd $input; ls *gz))
	myvar=0
	for id in ${sample[@]}
	do
		myvar=$(($myvar + 1 ))
		mkdir -p $output/1_qc/"${id%%.*}"
		### fastqc shell ###
		time fastqc -t 24 \
		-q $input/$id \
		--outdir $output/1_qc/"${id%%.*}" > $output/1_qc/"${id%%.*}"/${id%%.*}_qc.log 2>&1 &
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	wait
	echo '-----------------------------------------------------------'
	echo "end time: `date` | complete 1_qc"
	echo '-----------------------------------------------------------'
fi

if [ ! -d "$output/2_trim" ]; then
	echo '-----------------------------------------------------------'
	echo "start time: `date` | procress 2_trim"
	echo '-----------------------------------------------------------'
	mkdir -p $output/2_trim
	mkdir -p $output/2_trim/report
	myvar=0
	for id in ${sample_name[@]}
	do
		myvar=$(($myvar + 1))
		### fastp trim shell ###
		time fastp --thread 64 \
		--detect_adapter_for_pe \
		--overrepresentation_analysis \
		--trim_front1 2 --trim_front2 2 \
		--cut_front --cut_tail --cut_window_size 3 \
		--cut_mean_quality 30 \
		-i $input/${id}_1_trim.fastq -I $input/${id}_2_trim.fastq \
		-o $output/2_trim/${id}_1_trimed.fastq \
		-O $output/2_trim/${id}_2_trimed.fastq \
		-j $output/2_trim/report/${id}.json \
		-h $output/2_trim/report/${id}.html \
		-R "$id" >$output/2_trim/report/${id}_trim.log 2>&1 &
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	wait
	echo '--------------------------------------------------------'
	echo "end time: `date` | complete 2_trim"
	echo '--------------------------------------------------------'
fi

if [ ! -d "$output/3_align" ];then
	echo '-----------------------------------------------------------'
	echo "start time: `date` | procress 3_align"
	echo '-----------------------------------------------------------'
	mkdir -p $output/3_align
	myvar=0
	for id in ${sample_name[@]}
	do
		myvar=$(($myvar + 1))
		time STAR --runThreadN 48 \
		--genomeDir /opt/test/ref_index/star_index \
		--readFilesIn $output/2_trim/${id}_1_trimed.fastq \
		$output/2_trim/${id}_2_trimed.fastq \
		--outFileNamePrefix $output/3_align/${id}_ \
		--outSAMtype BAM Unsorted \
		--outSAMstrandField intronMotif \
		--outSAMattributes All \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outTmpDir $output/3_align/${id}_STAR_tmp \
		--twopassMode Basic \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 8 \
		--chimOutJunctionFormat 1 \
		--alignMatesGapMax 100000 \
		--alignIntronMax 100000 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 > $output/3_align/${id}_align.log 2>&1 &
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	wait
	echo '----------------------------------------------------------'
	echo " end time: `date` | complete 3_align "
	echo '----------------------------------------------------------'
fi

if [ ! -d "$output/4_sort" ];then
	echo '-----------------------------------------------------------'
	echo " start time: `date` | procress 4_sort "
	echo '-----------------------------------------------------------'
	mkdir -p $output/4_sort
	myvar=0
	for id in ${sample_name[@]}
	do 
		myvar=$(($myvar + 1))
		### BAM sort shell ###
		time samtools sort -@ 32 \
		-O bam \
		-o $output/4_sort/${id}_sorted.bam \
		$output/3_align/${id}_Aligned.out.bam
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	echo '-----------------------------------------------------------'
	echo " end time: `date` | complete 4_sort "
	echo '-----------------------------------------------------------'
fi

if [ ! -d "$output/5_qcbam" ];then
	echo '------------------------------------------------------------'
	echo " start time: `date` | procress 5_qcbam"
	echo '------------------------------------------------------------'
	mkdir -p $output/5_qcbam
	myvar=0
	for id in ${sample_name[@]}
	do
		myvar=$(($myvar + 1))
		### bam qullimap shell ###
		time qualimap rnaseq \
		-bam $output/4_sort/${id}_sorted.bam \
		-gtf /opt/test/xiongyifan/reference/hg38_chr.gtf \
		-pe -s\
		-outdir $output/5_qcbam/$id > $output/5_qcbam/${id}.bamqc.log 2>&1 
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	echo '----------------------------------------------------------'
	echo " end time: `date` | complete 5_qcbam"
	echo '----------------------------------------------------------'
fi

if [ ! -d "$output/6_counts" ];then
	echo '------------------------------------------------------------'
	echo "start time: `date` | procress 6_counts"
	echo '------------------------------------------------------------'
	mkdir -p $output/6_counts
	myvar=0
	for id in ${sample_name[@]}
	do
		myvar=$(($myvar + 1))
		### feature counts shell ###
		time featureCounts -T 64 \
		-Q 10 -F GTF\
		-t exon -g gene_id \
		--primary -s 0 -p \
		-a /opt/test/xiongyifan/reference/hg38_chr.gtf \
		-o $output/6_counts/${id}_counts \
		$output/4_sort/${id}_sorted.bam > $output/6_counts/${id}_counts.log 2>&1 
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	echo '------------------------------------------------------------'
	echo " end time: `date` | complete 6_counts"
	echo '------------------------------------------------------------'
fi

if [ ! -d "$output/7_fusion" ];then
	echo '-------------------------------------------------------------'
	echo "start time: `date` | procress 7_fusion"
	echo '-------------------------------------------------------------'
	mkdir -p $output/7_fusion
	myvar=0
	treatment=($(cd $input/; ls stroke*gz))
	lib_dir=/opt/test/xiongyifan/rnaseq/public_data_test/fusion/CTAT_resource_lib/ctat_genome_lib_build_dir
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time STAR-Fusion\
		--genome_lib_dir $lib_dir \
		-J $output/3_align/${id}_Chimeric.out.junction \
		--output_dir $output/7_fusion/${id} > $output/7_fusion/${id}_fusion.log 2>&1 &
		if [ "$myvar" = "6" ]
		then
			myvar=0
			wait
		fi
	done
	echo '-------------------------------------------------------------'
	echo "start time: `date` | complete 7_fusion"
	echo '-------------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel" ];then
	echo '--------------------------------------------------------------'
	echo "start time: `date` | procress 8_SNPindel"
	echo '--------------------------------------------------------------'
	mkdir -p $output/8_SNPindel
	mkdir -p $output/8_SNPindel/rmdup
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do
		time gatk --java-options "-Xmx6G"\
		MarkDuplicates --REMOVE_DUPLICATES true \
		--INPUT $output/4_sort/${id}_sorted.bam \
		--METRICS_FILE $output/8_SNPindel/rmdup/${id}_rmdup.txt \
		--OUTPUT $output/8_SNPindel/rmdup/${id}_rmdup.bam > $output/8_SNPindel/rmdup/${id}_rmdup.log 2>&1 &
		myvar=$(($myvar +1))
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
		wait
	done
	echo '--------------------------------------------------------'
	echo "end time: `date` | finshed rmdup"
	echo '--------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/splitN" ];then
	echo '--------------------------------------------------------'
	echo "end time: `date` | procress splitn"
	echo '--------------------------------------------------------'
	mkdir -p $output/8_SNPindel/splitN
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do 
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6G" \
		SplitNCigarReads \
		-R $hg38_seq \
		-I $output/8_SNPindel/rmdup/${id}_rmdup.bam \
		-O $output/8_SNPindel/splitN/${id}_splitn.bam >$output/8_SNPindel/splitN/${id}_splitn.log 2>&1 &
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
	done
	echo '-------------------------------------------------------'
	echo "end time: `date` | finshed splitn"
	echo '-------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/addRG" ];then
	echo '-------------------------------------------------------'
	echo "end time: `date` | procress addRG"
	echo '-------------------------------------------------------'
	mkdir -p $output/8_SNPindel/addRG
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6G" \
		AddOrReplaceReadGroups \
		-I $output/8_SNPindel/splitN/${id}_splitn.bam\
		-O $output/8_SNPindel/addRG/${id}_addRG.bam \
		-RGLB ${id} \
		-RGPL ILLUMINA \	
		-RGPU ${id} \
		-RGSM ${id}
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
		wait
	done
	echo '--------------------------------------------------------'
	echo "end time: `date` | finshed addRG"
	echo '--------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/recaltable" ];then
	echo '--------------------------------------------------------------'
	echo "end time: `date` | procress recaltable"
	echo '--------------------------------------------------------------'
	mkdir -p $output/8_SNPindel/recaltable
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6G"\
		BaseRecalibrator \
		-I $output/8_SNPindel/addRG/${id}_addRG.bam\
		-R $hg38_seq \
		--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/1000G_omni2.5.hg38.vcf \
		--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/1000G_phase1.snps.high_confidence.hg38.vcf  \
		--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/hapmap_3.3.hg38.vcf \
		--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/dbsnp_146.hg38.vcf \
		--known-sites /opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-O $output/8_SNPindel/recaltable/${id}_recal.table > $output/8_SNPindel/recaltable/${id}_table.log 2>&1 &
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
		wait
	done
	echo '--------------------------------------------------------------'
	echo "end time: `date` | finshed recaltable"
	echo '--------------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/recal" ];then
	echo '------------------------------------------------------'
	echo "end time: `date` | procress recal"
	echo '------------------------------------------------------'
	mkdir -p $output/8_SNPindel/recal
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6G"\
		ApplyBQSR \
		-I $output/8_SNPindel/addRG/${id}_addRG.bam \
		-R $hg38_seq \
		--bqsr-recal-file $output/8_SNPindel/recaltable/${id}_recal.table \
		-O $output/8_SNPindel/recal/${id}_recal.bam >$output/8_SNPindel/recal/${id}_recal.log 2>&1 &
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
		wait
	done
	echo '------------------------------------------------------'
	echo "end time: `date` | finshed recal"
	echo '------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/vcf" ];then
	echo '----------------------------------------------------------------'
	echo "end time: `date` | procress variant calling"
	echo '----------------------------------------------------------------'
	mkdir -p $output/8_SNPindel/vcf
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	hg38_seq=/opt/tsinghua/cfDNApipeTest/file/vcf/hg38/Homo_sapiens_assembly38.fasta
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6g" \
		HaplotypeCaller \
		-I $output/8_SNPindel/recal/${id}_recal.bam \
		-R $hg38_seq \
		-O $output/8_SNPindel/vcf/${id}.vcf \
		--dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
		wait
	done
	echo '----------------------------------------------------------------'
	echo "end time: `date` | finshed variant calling"
	echo '----------------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/1" ];then
	echo '----------------------------------------------------------'
	echo "end time: `date` | procress vcf index"
	echo '----------------------------------------------------------'
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk IndexFeatureFile \
		-I $output/8_SNPindel/vcf/${id}.vcf
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
	done
	echo '----------------------------------------------------------'
	echo "end time: `date` | finshed vcf index"
	echo '----------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/1" ];then
	echo '---------------------------------------------------------------'
	echo "end time: `date` | procress vcf filtration"
	echo '---------------------------------------------------------------'
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk --java-options "-Xmx6g" \
		VariantFiltration \
		-V $output/8_SNPindel/vcf/${id}.vcf \
		-R $hg38_seq \
		-O $output/8_SNPindel/vcf/${id}_filted.vcf\
		-window 35 \
		-cluster 3 \
		-filter "FS > 30.0 || QD < 2.0" \
		--filter-name FSQD 2
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
	done
	echo '---------------------------------------------------------------'
	echo "end time: `date` | finshed vcf filtration"
	echo '---------------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/1" ];then
	echo '-------------------------------------------------------------'
	echo "end time: `date` | procress filted index"
	echo '-------------------------------------------------------------'
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time gatk IndexFeatureFile \
		-I $output/8_SNPindel/vcf/${id}_filted.vcf
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
		fi
	done
	echo '-------------------------------------------------------------'
	echo "end time: `date` | finshed filted index"
	echo '-------------------------------------------------------------'
fi

if [ ! -d "$output/8_SNPindel/annovar" ];then
	echo '-----------------------------------------------------------'
	echo "end time: `date` | procress annovar"
	echo '-----------------------------------------------------------'
	mkdir -p $output/8_SNPindel/annovar
	myvar=0
	treatment=($(cd $input/; ls stroke*gz|sed 's/_..fastq.gz//g'|uniq))
	for id in ${treatment[@]}
	do
		myvar=$(($myvar +1))
		time /opt/test/pacbio/annovar/table_annovar.pl \
		$output/8_SNPindel/vcf/${id}_filted.vcf \
		/opt/test/pacbio/annovar/humandb \
		--out $output/8_SNPindel/annovar/${id}_annovar.vcf \
		-buildver hg38 -otherinfo -nastring . \
		-protocol refGeneName,refGene,mim2gene,Gencode,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp42a \
		--operation r,g,r,r,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f \
		--vcfinput --thread 30 --remove --intronhgvs 300 \
		--argument '--colsWanted 5,--transcript_function,,--colsWanted 5,,,,,,,,,,,,,,,,'
		if [ "myvar" = "3" ]
		then
			myvar=0
			wait
			echo '-----------------------------------------------------------'
			echo "end time: `date` | complete 8_SNPindel"
			echo '-----------------------------------------------------------'
		fi
	done
fi



















































































