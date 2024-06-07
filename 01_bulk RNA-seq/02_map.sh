#!/bin/bash

tar=/01_result/07_lung/03_bulkRNAseq/siFOXO3
raw=/01_rawdata/01_lung/bulkRNAseq/siFOXO3/Rawdata
map=$tar/03_map

##for sample in bulk-Ctrl-1 bulk-Ctrl-2 bulk-Ctrl-3 bulk-Ctrl-4 bulk-Ctrl-5 bulk-Ctrl-6 bulk-Ctrl-7 bulk-AL_TB-1 bulk-L_TB-1 bulk-AL_TB-2 bulk-L_TB-2 bulk-AL_TB-3 bulk-L_TB-3 bulk-Ctrl-8 bulk-Ctrl-9 bulk-Ctrl-10 bulk-Ctrl-11 bulk-Ctrl-12 bulk-Ctrl-13 bulk-AL_TB-4 bulk-AL_TB-5 bulk-AL_TB-6 bulk-AL_TB-7 bulk-AL_TB-8 bulk-AL_TB-9 bulk-AL_TB-10 bulk-AL_TB-11 bulk-AL_TB-12 bulk-AL_TB-13 bulk-AL_TB-14 bulk-AL_TB-15 bulk-AL_TB-16 bulk-AL_TB-17 bulk-AL_TB-18 bulk-L_TB-4 bulk-L_TB-5 bulk-L_TB-6 bulk-L_TB-7 bulk-L_TB-8 bulk-L_TB-9 bulk-L_TB-10 bulk-L_TB-11 bulk-L_TB-12 bulk-L_TB-13 bulk-L_TB-14 bulk-L_TB-15 bulk-L_TB-16 bulk-L_TB-17 bulk-L_TB-18 si-NC-1 si-NC-2 si-NC-3 si-FOXO3-1 si-FOXO3-2 si-FOXO3-3 Vehicle-1 Vehicle-2 Vehicle-3 Thrombin-1 Thrombin-2 Thrombin-3 WT_Vehicle-1 WT_Vehicle-2 WT_Vehicle-3 WT_Thrombin-1 WT_Thrombin-2 WT_Thrombin-3 p65KO_Thrombin-1 p65KO_Thrombin-2 p65KO_Thrombin-3;
for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$map/scripts/${sample}_map.sh
echo "#PBS -N map_$sample">>$map/scripts/${sample}_map.sh
echo "#PBS -o $sample.out">>$map/scripts/${sample}_map.sh
echo "#PBS -e $sample.err">>$map/scripts/${sample}_map.sh
echo "#PBS -l nodes=1:ppn=12">>$map/scripts/${sample}_map.sh
echo "#PBS -q batch">>$map/scripts/${sample}_map.sh
echo '#PBS -p 1023'>>$map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'
echo Mapping:$sample is Starting
trim=$tar/02_trim
map=$tar/03_map

clean1=$trim/$sample/*_R1_val_1.fq.gz
clean2=$trim/$sample/*_R2_val_2.fq.gz

index=/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/01_hisat2_index/hg19

result=$map/$sample
log=$map/logs
mkdir $result

hisat2=/anaconda3/bin/hisat2
$hisat2 -p 24 -x $index -1 $clean1 -2 $clean2 --dta -S $result/${sample}.sam 2>$log/${sample}_map.log

echo Mapping:$sample is Done

bam=$result/${sample}.bam
sort=$result/${sample}.sort.bam
sort1=$result/${sample}.sort.name.bam

samtools=/anaconda3/bin/samtools
$samtools view -bS -@ 24 -q 10 $result/${sample}.sam>$bam &&
$samtools sort -@ 24 $bam -o $sort &&
$samtools sort -@ 24 -n $bam -o $sort1 &&
$samtools index -@ 24 $sort


'>>$map/scripts/${sample}_map.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$map/run_map.sh

