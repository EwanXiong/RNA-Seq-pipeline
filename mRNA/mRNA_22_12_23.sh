##############################################################################
####                 RNA-seq standered analysis pipeline                  ####
####                          Author  EvanXiong                           ####
####                        Update time 2022.12.23                        ####  
##############################################################################
# Analysis including:
# 1. Rawdata quality control
# 2. Adapater trimming
# 3. Reference genome alignment 
# 4. Bam file sort
# 5. Bam file quality control
# 6. Read counts calling
# 7. Alternative splicing analysis
# 8. Novel transcripts identification

input=/opt/.../rawdata
output=/opt/.../output
index=/opt/.../star
gff=/opt/.../Genome.gff3
gtf=/opt/.../Genome.gtf
genome=/opt/.../Genome.fa

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
        -i $input/${id}_R1.fq.gz -I $input/${id}_R2.fq.gz \
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
        -outdir $output/5_qcbam/$id > $output/5_qcbam/${id}.log 2>&1 
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
    -Q 10 -F GFF \
    -t gene -g gene_id \
    --primary -s 0 -p \
    -a $gff \
    -o $output/6_counts/all_counts \
    $output/4_sort/*bam > $output/6_counts/featurecounts.log 2>&1 
    wait
    echo '------------------------------------------------------------'
    echo " end time: `date` | complete 6_counts"
    echo '------------------------------------------------------------'
fi

if [ ! -d "$output/7" ];then
    sed -i '/#/d' $output/6_counts/all_counts
    sed -i "s#$output/4_sort/##g" $output/6_counts/all_counts
    sed -i 's/_sorted.bam//g' $output/6_counts/all_counts
fi

if [ ! -d "$output/7_alternative" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress 7_alternative"
    echo '---------------------------------------------------------------'
    mkdir -p $output/7_alternative
    source ~/anaconda3/bin/activate rmats
    myvar=0
    cat group_id|while read id
    do
        myvar=$(($myvar + 1))
        time rmats.py \
        --b1 ${id%%_vs*}.txt \
        --b2 ${id#*vs_}.txt \
        --gtf $gtf \
        -t paired \
        --readLength 100 --variable-read-length \
        --cstat 0.0001 \
        --nthread 80 \
        --od $output/7_alternative/${id} \
        --tmp $output/7_alternative/${id}/tmp > $output/7_alternative/${id}/run.log 2>&1 &
        if [ "$myvar" = "6" ]
        then
            myvar=0
            wait
        fi
    done
    echo '---------------------------------------------------------------'
    echo "start time: `date` | complete 7_alternative"
    echo '---------------------------------------------------------------'
fi

if [ ! -d "$output/7_alternative/control" ]; then
    echo '---------------------------------------------------------------'
    echo "start time: `date` | procress alterna_plot"
    echo '---------------------------------------------------------------'
    mkdir -p $output/7_alternative/control
    source ~/anaconda3/bin/activate rmatsplot
    myvar=0
    cat group_id|while read i
    do 
        cat rmats_id |while read id
        do 
            time rmats2sashimiplot \
            --b1 `cat ${i#*vs_}.txt` \
            --b2 `cat ${i%%_vs*}.txt` \
            --l1 ${i#*vs_} \
            --l2 ${i%%_vs*} \
            -o $output/7_alternative/${i}/plot/plot_${id} \
            -t $id --exon_s 1 --intron_s 5 \
            -e $output/7_alternative/${i}/${id}_plot.txt &
        done
        wait
    done
    wait
    echo '---------------------------------------------------------------'
    echo "start time: `date` | complete alterna_plot"
    echo '---------------------------------------------------------------'
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
        -B -G $gff \
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
    -G $gff \
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
    -R -r $gff \
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
    source ~/anaconda3/bin/activate gffread
    mkdir -p $output/novel
    awk '{if ($3=="u"){print $0}}' $output/9_gffmerge/merged.stringtie_merged.gtf.tmap > $output/novel/filter1_by_u.tmap
    wait
    awk '($6>1 && $10>=200){print$0}' $output/novel/filter1_by_u.tmap > $output/novel/filter2_by_length.tmap
    wait
    awk '{print $5}' $output/novel/filter2_by_length.tmap > $output/novel/filter2_transcript_ID
    wait
    grep -w -Ff $output/novel/filter2_transcript_ID -w $output/9_gffmerge/merged.annotated.gtf > $output/novel/filter2_transcript.gtf
    wait
    awk '{if($3=="exon")print}' $output/novel/filter2_transcript.gtf > $output/novel/filter3_transcript_exon.gtf
    wait
    gffread -w $output/novel/filter3_transcript_exon.fa -g $genome $output/novel/filter3_transcript_exon.gtf
    wait
    /root/anaconda3/envs/cpc2/bin/python /opt/tsinghua/software/cpc2/CPC2-beta/bin/CPC2.py \
    -i $output/novel/filter3_transcript_exon.fa -o $output/novel/cpc2_result > $output/novel/cpc2.log 2>&1 &
    wait
    awk '{if($8=="coding")print $1}' $output/novel/cpc2_result.txt > $output/novel/novel_coding_id.txt
    wait
    cat $output/novel/filter3_transcript_exon.fa |/opt/tsinghua/software/seqkit/seqkit grep -f $output/novel/novel_coding_id.txt > $output/novel/novel_transcript.fa
    wait
    /opt/tsinghua/software/seqkit/seqkit translate $output/novel/novel_transcript.fa > $output/novel/novel_transcript_pep.fa
    wait
    grep -w -Ff $output/novel/novel_coding_id.txt -w $output/novel/filter2_transcript.gtf > $output/novel/novel_transcript.gtf
    wait
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
