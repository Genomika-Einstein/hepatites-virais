data="/home/ubuntu/bioinfo/data"

# recebe os inputs
while [ "$1" != "" ]
        do
        case $1 in
                "-s")
                        shift
                        sample="$1" # sample name (string)
                        shift
                        ;;
                "-1")
                        shift
                        read1="$1" # read 1 (.fastq)
                        shift
                        ;;
                "-2")
                        shift
                        read2="$1" # read 2 (.fastq)
                        shift
                        ;;
		"-ref")
			shift
			ref="$1"
			shift
			;;
        esac
        done

# alinhar reads
# OUTPUT FILES:
# [BWA] _bwamem.log
# [samtoools] .bam
echo "==============="
echo "Alinhando reads"
cur=$(pwd)
cd $(dirname $ref)
bwa mem -v 1 -M -R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" $(basename $ref).fasta $read1 $read2 \
2>$data/log/${sample}_bwamem.log | samtools view -F 4 -Sb ->$data/BAM/${sample}.bam 2>$data/log/${sample}_samtoolsview.log
cd $cur

# sort + index bam
# [samtools] .bam.bai
samtools sort $data/BAM/${sample}.bam $data/BAM/${sample}_sorted
mv $data/BAM/${sample}_sorted.bam $data/BAM/${sample}.bam
samtools index $data/BAM/${sample}.bam $data/BAM/${sample}.bam.bai
