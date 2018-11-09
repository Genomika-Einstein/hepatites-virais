afterqc="/home/ubuntu/bioinfo/app/AfterQC-0.9.7"
fastqc="/home/ubuntu/bioinfo/app/FastQC/"

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
        esac
        done

if [ ! -d ../data/FastQC/${sample}_fastqc ]
	then mkdir ../data/FastQC/${sample}_fastqc
fi
echo Iniciando FastQC R1
$fastqc/fastqc -o ../data/FastQC/${sample}_fastqc -q $read1
echo Iniciando FastQC R2
$fastqc/fastqc -o ../data/FastQC/${sample}_fastqc -q $read2
echo Iniciando AfterQC
afterqc_report="/home/ubuntu/bioinfo/data/AfterQC/"$sample
filtered_reads="/home/ubuntu/bioinfo/data/fastq_filtered/"$sample
pypy $afterqc/after.py -1 $read1 -2 $read2 -r $afterqc_report --trim_front=-1 --trim_tail=-1 --seq_len_req=100 -g $filtered_reads -b .
rm ${sample}*.bad.fq.gz
