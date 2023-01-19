# Trim data extract for report writting
for i in `ls *html`;do echo ${i%%.*} >> sample.txt;done

cat *log|grep -A 4 'after filtering'|grep 'total reads:'|sed 's/total reads: //g' >R1.txt

cat R1.txt |awk '{print $1 * 2 / 1000000}' > Total_clean_reads.txt

cat *log|grep -A 4 'after filtering'|grep 'total bases:'|sed 's/total bases: //g' > R1_base.txt

cat *log|grep -A 4 'aftering filtering'|grep 'total bases:'|sed 's/total bases: //g' > R2_base.txt

paste R1_base.txt R2_base.txt > all_base.txt

cat all_base.txt |awk '{print ($1 + $2) / 1000000000 }' > Total_base.txt

cat *json | grep -w -A 9 'after_filtering' |grep 'gc_content' |sed 's/.*://g' |awk '{print $1 *100}' > GC.txt

for i in `cat *log|grep -A 4 'aftering filtering'|grep 'Q20 bases:'|sed 's/Q20 bases: //g'`;do echo ${i#*(} |sed 's/%)//g' >> Q20.txt;done 

for i in `cat *log|grep -A 4 'aftering filtering'|grep 'Q30 bases:'|sed 's/Q30 bases: //g'`;do echo ${i#*(} |sed 's/%)//g' >> Q30.txt;done 

paste sample.txt Total_clean_reads.txt Total_base.txt GC.txt Q20.txt Q30.txt > trim_sum.txt
















