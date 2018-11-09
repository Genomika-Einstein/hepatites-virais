sample_list="$(cat samples.txt)"
#sample_list="9648_Hepatite_S2 AAS_IMT_HEPATITE_S7 RPI63_IMT_HEPATITE_S3"
data="../data"
ref="/home/ubuntu/bioinfo/reference/HBV_ref/HBV_ref"
nuc_res="/home/ubuntu/bioinfo/reference/HBV_RT.txt"

for sample in $sample_list
do
	echo Iniciando amostra $sample
	# Encontra read1 e read2
	cur=$(pwd)
	cd $data/fastq
	var=$(find `pwd` -type f | grep $sample | grep .fastq)
	cd $cur
	read1=$(echo $var | cut -d' ' -f1)
	read2=$(echo $var | cut -d' ' -f2)
	if [ $read1 \> $read2 ]
	then
		aux=$read1
		read1=$read2
		read2=$aux
	fi
	# controle de qualidade
	#bash quality_control.sh -s $sample -1 $read1 -2 $read2
	read1=$(bash abs.sh $data/fastq_filtered/$sample/${sample}_R1_001.good.fq.gz)
	read2=$(bash abs.sh $data/fastq_filtered/$sample/${sample}_R2_001.good.fq.gz)
	# Assembly + Genotipagem
	bash assembly.sh -s $sample -1 $read1 -2 $read2
	# checagem de contaminação humana
	#bash human_contamination.sh $sample
	# alinhar contra referencia
	#bash bwa_call.sh -s $sample -1 $read1 -2 $read2 -ref $ref
	# chamar variantes
	#bash var_call.sh -s $sample -ref $ref
	# anotar variantes
	#python3 vcf_annotator.py $data/VCF/$sample.vcf ${ref}_features.fasta $ref.bed $data/variants/$sample
	#grep -Fx -f $data/variants/$sample.ptm.txt $nuc_res > $data/variants/$sample.nuc.txt
	# visualização no IGV
	#bash igv_call.sh -sample $sample
done
