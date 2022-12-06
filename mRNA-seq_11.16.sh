input=/opt/tsinghua/NuoHe/mRNA_22.11.17/rawdata
output=/opt/tsinghua/NuoHe/mRNA_22.11.17/output
index=/opt/tsinghua/cfDNApipeTest/file/mm/mm10/star
gtf=/opt/tsinghua/cfDNApipeTest/file/mm/mm10/mm10_chr.gtf
b1=/opt/tsinghua/NuoHe/mRNA_22.11.17/shell/ko.txt
b2=/opt/tsinghua/NuoHe/mRNA_22.11.17/shell/wt.txt
genome=/opt/tsinghua/cfDNApipeTest/file/mm/mm10/mm10.fa
# 7_Alternative use bam path 
bam_path_1=`cat ko.txt`
bam_path_2=`cat wt.txt`
# Alternative plot use group 
group_1=WT
group_2=KO

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
        time /opt/tsinghua/software/fastqc/fastqc/FastQC/fastqc -t 80 \
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
    cat id |while read id 
    do
        myvar=$(($myvar + 1))
        time fastp --thread 60 \
        --detect_adapter_for_pe \
        --overrepresentation_analysis \
        --trim_front1 2 --trim_front2 2 \
        --cut_front --cut_tail --cut_window_size 3 \
        --cut_mean_quality 30 \
        -i $input/${id}_1.fq.gz -I $input/${id}_2.fq.gz \
        -o $output/2_trim/${id}_1_trimed.fq.gz \
        -O $output/2_trim/${id}_2_trimed.fq.gz \
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
    cat id |while read id 
    do
        myvar=$(($myvar + 1))
        time STAR --runThreadN 60 \
        --genomeDir $index \
        --readFilesCommand zcat \
        --readFilesIn $output/2_trim/${id}_1_trimed.fq.gz \
        $output/2_trim/${id}_2_trimed.fq.gz \
        --outFileNamePrefix $output/3_align/${id}_ \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMattributes All \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outTmpDir $output/3_align/${id}_STAR_tmp > $output/3_align/${id}_align.log 2>&1 &
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
    cat id |while read id
    do 
        myvar=$(($myvar + 1))
        time samtools sort -@ 60 \
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
    cat id |while read id
    do
        myvar=$(($myvar + 1))
        time qualimap rnaseq \
        -bam $output/4_sort/${id}_sorted.bam \
        -gtf $gtf \
        -pe -s \
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
    time featureCounts -T 64 \
    -Q 10 -F GTF \
    -t exon -g gene_id \
    --primary -s 0 -p \
    -a $gtf \
    -o $output/6_counts/all_counts \
    $output/4_sort/*bam > $output/6_counts/featurecounts.log 2>&1 
    wait
    echo '------------------------------------------------------------'
    echo " end time: `date` | complete 6_counts"
    echo '------------------------------------------------------------'
fi

if [ ! -d "$output/7" ];then
    sed -i '/#/d' $output/6_counts/all_counts
    sed -i 's#/opt/tsinghua/NuoHe/mRNA_22.11.14/output/4_sort/##g' $output/6_counts/all_counts
    sed -i 's/_sorted.bam//g' $output/6_counts/all_counts
fi

if [ ! -d "$output/7_Alternative" ]; then
    echo '-----------------------------------------------------------'
    echo "start time: `date` | procress 7_alter"
    echo '-----------------------------------------------------------'
    mkdir -p $output/7_Alternative
    source ~/anaconda3/bin/activate rmats
    wait
    time rmats.py --b1 $b1 --b2 $b2 \
    --gtf $gtf \
    -t paired --readLength 100 --variable-read-length \
    --cstat 0.0001 \
    --nthread 80 \
    --od $output/7_Alternative \
    --tmp $output/7_Alternative > $output/7_Alternative/run.log &
    wait
    echo '-----------------------------------------------------------'
    echo "start time: `date` | complete 7_alter"
    echo '-----------------------------------------------------------'
fi

cat rmats_id |while read id
do head -n 11 $output/7_Alternative/${id}.MATS.JC.txt > $output/7_Alternative/${id}_plot.txt
done

if [ ! -d "$output/7_Alternative/plot" ]; then
     echo '-----------------------------------------------------------'
     echo "start time: `date` | procress plot_alter"
     echo '-----------------------------------------------------------'
    mkdir -p $output/7_Alternative/plot
    source ~/anaconda3/bin/activate rmatsplot
    cat rmats_id |while read id
    do rmats2sashimiplot --b1 $bam_path_1 --b2 $bam_path_2 \
    --l1 $group_1 --l2 $group_2 \
    -o $output/7_Alternative/plot/plot_${id} \
    -t ${id} --exon_s 1 --intron_s 5 \
    -e $output/7_Alternative/${id}_plot.txt > $output/7_Alternative/plot/${id}.log 2>&1 &
    wait
    done
    echo '-----------------------------------------------------------'
    echo "start time: `date` | complete plot_alter"
    echo '-----------------------------------------------------------'
fi

if [ ! -d "$output/8_stringtie" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress 8_stringtie"
    echo '---------------------------------------------------------------'
    mkdir -p $output/8_stringtie
    source ~/anaconda3/bin/activate gffread
    myvar=0
    cat id |while read id
    do
        myvar=$(($myvar + 1))
        time stringtie -p 64 \
        -B -G $gtf \
        $output/4_sort/${id}_sorted.bam \
        -o $output/8_stringtie/${id}.stringtie.gtf \
        -A $output/8_stringtie/${id}.tab \
        -l ${id} > $output/8_stringtie/${id}_assmebly.log 2>&1 &
        if [ "$myvar" = "6" ]
        then
            myvar=0
            wait
        fi
    done
    wait
    echo '---------------------------------------------------------------'
    echo "end time: `date` | complete 8_stringtie"
    echo '---------------------------------------------------------------'
fi

if [ ! -d "$output/9_gffmerge" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress 9_gffmerge"
    echo '---------------------------------------------------------------'
    mkdir -p $output/9_gffmerge
    ls $output/8_stringtie/*gtf > $output/8_stringtie/mergelist.txt
    time stringtie --merge -p 64 \
    -G $gtf \
    -o $output/9_gffmerge/stringtie_merged.gtf \
    $output/8_stringtie/mergelist.txt > $output/9_gffmerge/merge.log 2>&1 & 
    wait
    echo '---------------------------------------------------------------'
    echo "end time: `date` | complete 9_gffmerge"
    echo '---------------------------------------------------------------'
    wait
fi

if [ ! -d "$output/10" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress gffcompare"
    echo '---------------------------------------------------------------'
    time gffcompare \
    -R -r $gtf \
    -o $output/9_gffmerge/merged \
    $output/9_gffmerge/stringtie_merged.gtf > $output/9_gffmerge/gffcompare.log 2>&1 &
    wait
    echo '---------------------------------------------------------------'
    echo "end time: `date` | complete gffcompare"
    echo '---------------------------------------------------------------'
fi

if [ ! -d "$output/novel" ]; then
    echo '--------------------------------------------------------------------'
    echo "start time: `date` | procress novel_transcript"
    echo '--------------------------------------------------------------------'
    mkdir -p $output/novel
    awk '{if ($3=="u"){print $0}}' $output/9_gffmerge/merged.stringtie_merged.gtf.tmap > $output/novel/filter1_by_u.tmap
    awk '($6>1 && $10>=200){print$0}' $output/novel/filter1_by_u.tmap > $output/novel/filter2_by_length.tmap
    awk '{print $5}' $output/novel/filter2_by_length.tmap > $output/novel/filter2_transcript_ID
    grep -w -Ff $output/novel/filter2_transcript_ID -w $output/9_gffmerge/merged.annotated.gtf > $output/novel/filter2_transcript.gtf
    awk '{if($3=="exon")print}' $output/novel/filter2_transcript.gtf > $output/novel/filter3_transcript_exon.gtf
    gffread -w $output/novel/filter3_transcript_exon.fa -g $genome $output/novel/filter3_transcript_exon.gtf
    /root/anaconda3/envs/cpc2/bin/python /opt/tsinghua/software/cpc2/CPC2-beta/bin/CPC2.py \
    -i $output/novel/filter3_transcript_exon.fa -o $output/novel/cpc2_result > $output/novel/cpc2.log 2>&1 &
    awk '{if($8=="coding")print $1}' $output/novel/cpc2_result.txt > $output/novel/novel_coding_id.txt
    cat $output/novel/filter3_transcript_exon.fa |/opt/tsinghua/software/seqkit/seqkit grep -f $output/novel/novel_coding_id.txt > $output/novel/novel_transcript.fa
    /opt/tsinghua/software/seqkit/seqkit translate $output/novel/novel_transcript.fa > $output/novel/novel_transcript_pep.fa
    grep -w -Ff $output/novel/novel_coding_id.txt -w $output/novel/filter2_transcript.gtf > $output/novel/novel_transcript.gtf
    emapper.py --cpu 60 -m diamond \
    --data_dir /opt/tsinghua/software/eggdb \
    --annotate_hits_table emapper_anno --no_file_comments \
    --dmnd_db /opt/tsinghua/software/eggdb/eggnog_proteins.dmnd \
    -i $output/novel/novel_transcript_pep.fa \
    -o $output/novel/novel --output_dir $output/novel/ --override &
    wait
    echo '---------------------------------------------------------------------'
    echo "end time: `date` | complete novel_transcript"
    echo '---------------------------------------------------------------------'
fi



去掉featurecount结果中的下游分析无关列
Rscript /opt/test/xiongyifan/script_R/trim_count_matrix.r \
--outputdir /opt/tsinghua/NuoHe/mRNA_22.11.24/output/6_counts \
--input all_counts

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
--config /opt/tsinghua/test_file/transform/Finalcd /config.yml \
--KEGG 
####GSEA用法
nohup Rscript /opt/tsinghua/test_file/transform/Final/GSEA.r \
--config /opt/tsinghua/test_file/transform/Final/config.yml \
--GESA


Rscript /opt/tsinghua/test_file/transform/Final/convert_ID.r --outputdir ./ \
--input /opt/tsinghua/NuoHe/mRNA_22.11.17/output/DEGs/all_counts.txt \
--type mm 



output=/opt/tsinghua/NuoHe/mRNA_22.11.24/output

if [ ! -d "$output/DEGs/output" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress DEGs"
    echo '---------------------------------------------------------------'
    myvar=0
    cat id |while read id
    do
        mkdir -p $output/DEGs/$id
        myvar=$(($myvar + 1))
        Rscript /opt/tsinghua/test_file/transform/Final/Auto_Enrichment_Final.r \
        --config ${id}.yml \
        --ALL > $output/DEGs/$id/run.log 2>&1 &
        if [ "$myvar" = "6" ]
        then
            myvar=0
            wait
        fi
    done
    wait
    echo '---------------------------------------------------------------'
    echo "end time: `date` | complete DEGs"
    echo '---------------------------------------------------------------'
fi











source ~/anaconda3/bin/activate rmats

    
rmats.py --b1 $b1 --b2 $b2 \
--gtf $gtf \
-t paired --readLength 100 --variable-read-length \
--cstat 0.0001 \
--nthread 80 \
--od $output/7_Alternative \
--tmp $output/7_Alternative > $output/7_Alternative/run.log &





if [ ! -d "$output/7_" ]; then
    echo '-----------------------------------------------------------'
    echo "start time: `date` | procress 7_alter"
    echo '-----------------------------------------------------------'
    source ~/anaconda3/bin/activate rmats
    wait
    cat bam_path |while read i
    do mkdir 

        time rmats.py ${i} \
    --gtf $gtf \
    -t paired --readLength 100 --variable-read-length \
    --cstat 0.0001 \
    --nthread 80 \
    --od $output/7_Alternative \
    --tmp $output/7_Alternative > $output/7_Alternative/run.log &
    wait
    echo '-----------------------------------------------------------'
    echo "start time: `date` | complete 7_alter"
    echo '-----------------------------------------------------------'
fi

rmats.py \
--b1 NoUV_Osmoter_L.txt --b2 NoUV_Osmoter_H.txt \
--gtf /opt/tsinghua/cfDNApipeTest/file/hg19/hg19.refGene.gtf \
-t paired --readLength 100 --variable-read-length \
--cstat 0.0001 \
--nthread 80 \
--od ./NoUV_Osmoter_H_vs_NoUVOsmoter_L \
--tmp ./NoUV_Osmoter_H_vs_NoUVOsmoter_L/tmp > ./NoUV_Osmoter_H_vs_NoUVOsmoter_L/run.log &


cat rmats_id |while read id
do head -n 11 $output/7_Alternative/${id}.MATS.JC.txt > $output/7_Alternative/${id}_plot.txt
done

if [ ! -d "$output/7_Alternative/plot" ]; then
     echo '-----------------------------------------------------------'
     echo "start time: `date` | procress plot_alter"
     echo '-----------------------------------------------------------'
    mkdir -p $output/7_Alternative/plot
    source ~/anaconda3/bin/activate rmatsplot
    cat rmats_id |while read id
    do rmats2sashimiplot --b1 $bam_path_1 --b2 $bam_path_2 \
    --l1 $group_1 --l2 $group_2 \
    -o $output/7_Alternative/plot/plot_${id} \
    -t ${id} --exon_s 1 --intron_s 5 \
    -e $output/7_Alternative/${id}_plot.txt > $output/7_Alternative/plot/${id}.log 2>&1 &
    wait
    done
    echo '-----------------------------------------------------------'
    echo "start time: `date` | complete plot_alter"
    echo '-----------------------------------------------------------'
fi


cat rmats_id|while read id 
do cat id |while read i 
do head -n 11 ${i}/${id}.MATS.JC.txt > ${i}/${id}_plot.txt;done
done


cat id |while read i 
do mkdir ${i}/plot; 


cat rmats_id|while read id 
do rmats2sashimiplot --b1 `cat NoUV_Osmoter_L.txt` --b2 `cat NoUV_Osmoter_H.txt` \
    --l1 NoUV_Osmoter_L --l2 NoUV_Osmoter_H \
    -o ./NoUV_Osmoter_H_vs_NoUVOsmoter_L/plot/plot_${id} \
    -t ${id} --exon_s 1 --intron_s 5 \
    -e ./NoUV_Osmoter_H_vs_NoUVOsmoter_L/${id}._plot.txt > ./NoUV_Osmoter_H_vs_NoUVOsmoter_L/plot/${id}.log 2>&1 &
done




报告信息统计

1. clean reads的质量
Read1 after filtering total reads:
cat *log|grep -A 4 'after filtering'|grep 'total reads:'|sed 's/total reads: //g'

Read2 after filtering:
cat *log|grep -A 4 'aftering filtering'|grep 'total reads:'|sed 's/total reads: //g'

Read1 after filtering total bases
cat *log|grep -A 4 'after filtering'|grep 'total bases:'|sed 's/total bases: //g'

Read2 after filtering total bases
cat *log|grep -A 4 'aftering filtering'|grep 'total bases:'|sed 's/total bases: //g'

Q20 bases: 
cat *log|grep -A 4 'aftering filtering'|grep 'Q20 bases:'|sed 's/Q20 bases: //g'

Q30 bases:
