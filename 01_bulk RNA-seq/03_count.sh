#!/bin/bash
tar=/01_result/07_lung/03_bulkRNAseq/siFOXO3
raw=/01_rawdata/01_lung/bulkRNAseq/siFOXO3/Rawdata
count=$tar/04_count

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$count/scripts/${sample}_count.sh
echo "#PBS -N count_$sample">>$count/scripts/${sample}_count.sh
echo "#PBS -o $sample.out">>$count/scripts/${sample}_count.sh
echo "#PBS -e $sample.err">>$count/scripts/${sample}_count.sh
echo "#PBS -l nodes=1:ppn=12">>$count/scripts/${sample}_count.sh
echo "#PBS -q batch">>$count/scripts/${sample}_count.sh
echo '#PBS -p 1023'>>$count/scripts/${sample}_count.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

gtf=/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87.final.gtf

uni=$tar/03_map
srt_bam=$uni/$sample/${sample}.sort.bam
count=$tar/04_count
result=$count/$sample
log=$count/logs
mkdir $result

echo HTseq:$sample is Starting

htseq=/anaconda3/bin/htseq-count
$htseq -f bam -r name -s no -a 10 $srt_bam $gtf > $result/${sample}_all.txt 2>$log/${sample}_count.log

echo HTseq-count has been Done'>>$count/scripts/${sample}_count.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$count/run_count.sh

