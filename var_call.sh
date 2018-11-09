data="../data" # path PRECISA ser relativo p/ HC
vt="../app/vt/vt"
while [ "$1" != "" ]
        do
        case $1 in
                "-s")
                        shift
                        sample="$1" # sample name (string)
                        shift
                        ;;
		"-ref")
			shift
			ref="$1" # ref path (fasta)
			shift
			;;
	esac
	done
#echo "Iniciando HaplotypeCaller"
#for ploidy in 1 5 10 20
#	do time gatk HaplotypeCaller -R $ref.fasta -I $data/BAM/$sample.bam -ploidy $ploidy -O $data/VCF/$sample.HC$ploidy.vcf
#done
#echo "Iniciando MUTECT2"
#time gatk Mutect2 -R $ref.fasta -I $data/BAM/$sample.bam -O $data/VCF/$sample.M2.vcf -tumor $sample
#$vt decompose -s $data/VCF/$sample.M2.vcf -o $data/VCF/$sample.vt.vcf
#mv $data/VCF/$sample.vt.vcf $data/VCF/$sample.M2.vcf
echo "Iniciando freebayes"
time freebayes -f $ref.fasta --pooled-continuous -p 1 -C 20 $data/BAM/$sample.bam | vcffilter -f "QUAL > 20" > $data/VCF/$sample.vcf
$vt decompose -s $data/VCF/$sample.vcf -o $data/VCF/$sample.vt.vcf
mv $data/VCF/$sample.vt.vcf $data/VCF/$sample.vcf
