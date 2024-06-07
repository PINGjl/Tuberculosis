#!/bin/bash
tar=/01_result/07_lung/03_bulkRNAseq/siFOXO3
raw=/01_rawdata/01_lung/bulkRNAseq/siFOXO3/Rawdata
count=$tar/04_count

##for sample in bulk-Ctrl-1 bulk-Ctrl-2 bulk-Ctrl-3 bulk-Ctrl-4 bulk-Ctrl-5 bulk-Ctrl-6 bulk-Ctrl-7 bulk-AL_TB-1 bulk-L_TB-1 bulk-AL_TB-2 bulk-L_TB-2 bulk-AL_TB-3 bulk-L_TB-3 bulk-Ctrl-8 bulk-Ctrl-9 bulk-Ctrl-10 bulk-Ctrl-11 bulk-Ctrl-12 bulk-Ctrl-13 bulk-AL_TB-4 bulk-AL_TB-5 bulk-AL_TB-6 bulk-AL_TB-7 bulk-AL_TB-8 bulk-AL_TB-9 bulk-AL_TB-10 bulk-AL_TB-11 bulk-AL_TB-12 bulk-AL_TB-13 bulk-AL_TB-14 bulk-AL_TB-15 bulk-AL_TB-16 bulk-AL_TB-17 bulk-AL_TB-18 bulk-L_TB-4 bulk-L_TB-5 bulk-L_TB-6 bulk-L_TB-7 bulk-L_TB-8 bulk-L_TB-9 bulk-L_TB-10 bulk-L_TB-11 bulk-L_TB-12 bulk-L_TB-13 bulk-L_TB-14 bulk-L_TB-15 bulk-L_TB-16 bulk-L_TB-17 bulk-L_TB-18 si-NC-1 si-NC-2 si-NC-3 si-FOXO3-1 si-FOXO3-2 si-FOXO3-3 Vehicle-1 Vehicle-2 Vehicle-3 Thrombin-1 Thrombin-2 Thrombin-3 WT_Vehicle-1 WT_Vehicle-2 WT_Vehicle-3 WT_Thrombin-1 WT_Thrombin-2 WT_Thrombin-3 p65KO_Thrombin-1 p65KO_Thrombin-2 p65KO_Thrombin-3;
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

