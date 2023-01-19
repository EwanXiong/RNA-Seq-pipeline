# alignment result sum for report writting
# need to prepare sampel.txt according to file name

cat *txt |grep 'Average mapped length'|sed 's/Average mapped length |//g' | sed 's/ //g'| sed 's/\t//g' > length.txt

cat *txt |grep 'Uniquely mapped reads number'|sed 's/Uniquely mapped reads number |//g'| sed 's/ //g'| sed 's/\t//g' > Unique_number.txt

cat *txt |grep 'Uniquely mapped reads %'|sed 's/Uniquely mapped reads % |//g' |sed 's/%//g'| sed 's/ //g'| sed 's/\t//g' > Unique_map.txt

cat *txt |grep 'Number of reads mapped to multiple loci'|sed 's/Number of reads mapped to multiple loci |//g'| sed 's/ //g'| sed 's/\t//g' >multi_reads.txt

cat *txt |grep '% of reads mapped to multiple loci'|sed 's/% of reads mapped to multiple loci |//g'|sed 's/%//g'| sed 's/ //g'| sed 's/\t//g' > multi_map.txt

cat *txt |grep 'Number of reads unmapped: too short'|sed 's/Number of reads unmapped: too short |//g' | sed 's/ //g'| sed 's/\t//g' > unmap_reads.txt

cat *txt |grep '% of reads unmapped: too short'|sed 's/% of reads unmapped: too short |//g'| sed 's/%//g'| sed 's/ //g'| sed 's/\t//g' >unmap.txt

paste sample.txt length.txt Unique_number.txt Unique_map.txt multi_reads.txt multi_map.txt unmap_reads.txt unmap.txt > sum_align.txt

