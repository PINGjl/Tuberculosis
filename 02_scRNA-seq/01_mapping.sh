count=/01_result/07_lung/01_ref

cd $count
for sample in Ctrl-1 Ctrl-2 Ctrl-3 Ctrl-4 Ctrl-5 Ctrl-6 Ctrl-7 AL-TB-1 AL-TB-2 AL-TB-3 L-TB-1 L-TB-2 L-TB-3;
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
