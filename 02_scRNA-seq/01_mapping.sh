
###Keep Mean reads closed 
cd /01_rawdata/01_lung_TB_new/data/L_TB-11
/02_Software/07_seqtk/seqtk/seqtk sample -2 -s100 TB_CD_B_17_S4_L005_R1_001.fastq.gz 0.6 > TB_CD_B_17_S60_L005_R1_001.fastq.gz
/02_Software/07_seqtk/seqtk/seqtk sample -2 -s100 TB_CD_B_17_S4_L005_R2_001.fastq.gz 0.6 > TB_CD_B_17_S60_L005_R2_001.fastq.gz

cd /dellstorage06/quj_lab/pingjiale/01_rawdata/01_lung_TB_new/data/AL_TB-11
/02_Software/07_seqtk/seqtk/seqtk sample -2 -s100 TB_CD_N_17_S2_L001_R1_001.fastq.gz 0.5 > TB_CD_N_17_S50_L005_R1_001.fastq.gz
/02_Software/07_seqtk/seqtk/seqtk sample -2 -s100 TB_CD_N_17_S2_L001_R2_001.fastq.gz 0.5 > TB_CD_N_17_S50_L005_R2_001.fastq.gz


###mapping
count=/01_result/07_lung/01_ref

cd $count
for sample in Ctrl-1 Ctrl-3 Ctrl-2 Ctrl-4 Ctrl-5 Ctrl-6 Ctrl-7 AL_TB-1 L_TB-1 AL_TB-2 L_TB-2 AL_TB-3 L_TB-3 Ctrl-8 Ctrl-9 Ctrl-10 Ctrl-11 Ctrl-12 Ctrl-13 AL_TB-4 AL_TB-5 AL_TB-6 AL_TB-7 AL_TB-8 AL_TB-9 AL_TB-10 AL_TB-11 AL_TB-12 L_TB-4 L_TB-5 L_TB-6 L_TB-7 L_TB-8 L_TB-9 L_TB-10 L_TB-11 L_TB-12;
do
/02_Software/cellranger/cellranger-4.0.0/bin/cellranger count \
--id=$sample \
--transcriptome=/03_Database/01_cellranger_4/tuberculosis/h19_tuberculosis/hg19_mtb \
--fastqs=/01_rawdata/05_Lung/01_data/$sample \
--sample=$sample \
--localcores=28 \
--localmem=100 \
--mempercore=300
#rm $count/$sample/outs/*.bam
done
