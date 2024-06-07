#!/bin/bash

tar=/01_result/07_lung/03_bulkRNAseq
raw=/01_rawdata/01_lung/bulkRNAseq/Rawdata
trim=$tar/02_trim
##for sample in bulk-Ctrl-1 bulk-Ctrl-2 bulk-Ctrl-3 bulk-Ctrl-4 bulk-Ctrl-5 bulk-Ctrl-6 bulk-Ctrl-7 bulk-AL_TB-1 bulk-L_TB-1 bulk-AL_TB-2 bulk-L_TB-2 bulk-AL_TB-3 bulk-L_TB-3 bulk-Ctrl-8 bulk-Ctrl-9 bulk-Ctrl-10 bulk-Ctrl-11 bulk-Ctrl-12 bulk-Ctrl-13 bulk-AL_TB-4 bulk-AL_TB-5 bulk-AL_TB-6 bulk-AL_TB-7 bulk-AL_TB-8 bulk-AL_TB-9 bulk-AL_TB-10 bulk-AL_TB-11 bulk-AL_TB-12 bulk-AL_TB-13 bulk-AL_TB-14 bulk-AL_TB-15 bulk-AL_TB-16 bulk-AL_TB-17 bulk-AL_TB-18 bulk-L_TB-4 bulk-L_TB-5 bulk-L_TB-6 bulk-L_TB-7 bulk-L_TB-8 bulk-L_TB-9 bulk-L_TB-10 bulk-L_TB-11 bulk-L_TB-12 bulk-L_TB-13 bulk-L_TB-14 bulk-L_TB-15 bulk-L_TB-16 bulk-L_TB-17 bulk-L_TB-18 si-NC-1 si-NC-2 si-NC-3 si-FOXO3-1 si-FOXO3-2 si-FOXO3-3 Vehicle-1 Vehicle-2 Vehicle-3 Thrombin-1 Thrombin-2 Thrombin-3 WT_Vehicle-1 WT_Vehicle-2 WT_Vehicle-3 WT_Thrombin-1 WT_Thrombin-2 WT_Thrombin-3 p65KO_Thrombin-1 p65KO_Thrombin-2 p65KO_Thrombin-3;
for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$trim/scripts/${sample}_trim.sh
echo "#PBS -N Trim_$sample">>$trim/scripts/${sample}_trim.sh
echo "#PBS -o $sample.out">>$trim/scripts/${sample}_trim.sh
echo "#PBS -e $sample.err">>$trim/scripts/${sample}_trim.sh
echo "#PBS -l nodes=1:ppn=28">>$trim/scripts/${sample}_trim.sh
echo "#PBS -q batch">>$trim/scripts/${sample}_trim.sh
echo '#PBS -p 1023'>>$trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming:$sample is Starting

trim=$tar/02_trim
log=$trim/logs

fq1=$raw/$sample/*_R1.fq.gz
fq2=$raw/$sample/*_R2.fq.gz

result=$trim/$sample
mkdir $result

trim_galore=/anaconda3/envs/py2/bin/trim_galore
$trim_galore --fastqc --path_to_cutadapt /anaconda3/envs/py2/bin/cutadapt --stringency 3 --paired --output_dir $result $fq1 $fq2 2>>$log/${sample}.log

echo Trimming has been Done'>>$trim/scripts/${sample}_trim.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$trim/run_trim.sh
