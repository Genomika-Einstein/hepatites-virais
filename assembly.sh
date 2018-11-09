# SINTAXE: bash assembly.sh -s sample -1 read1 -2 read2
# 1. Faz o assembly (spades)
# 2. Filtra o Scaffold
# 3. Blasta contra genomas de referencia
# 4. Genotipagem e assembly

spades="/home/ubuntu/bioinfo/app/SPAdes-3.11.1-Linux/bin/spades.py" # path para o spades
data="/home/ubuntu/bioinfo/data"
reference="/home/ubuntu/bioinfo/reference/hbv_all_genotypes/hbv_all_genotypes.fasta" # path para os genomas de referencia

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

# 1. Assembly
#if [ ! -d ${sample}_spades ]
#then
#	mkdir ${sample}_spades
#fi

#python $spades -1 $read1 -2 $read2 -o $data/spades/${sample}_spades -k 21,33,55,77,99,127 --phred-offset 33 -m 1 -t 1 --only-assembler
##python $spades -1 $read1 -2 $read2 -o $data/spades/${sample}_spades -k 21 --phred-offset 33 -m 1 --only-assembler -t 1
python3 filter_scaffold.py $data/spades/${sample}_spades/scaffolds.fasta $data/contigs/${sample}_contigs.fasta
blastn -query $data/contigs/${sample}_contigs.fasta -subject $reference -word_size 10 -outfmt '6 std sstrand' -max_target_seqs 3 > $data/blastn/${sample}_blast.txt
cp $data/blastn/${sample}_blast.txt $data/blastn/${sample}_blast.filtered.txt
python3 process_blast.py $data/blastn/${sample}_blast.filtered.txt $data/contigs/${sample}_contigs.fasta $sample > $data/WGS/${sample}_scaffold.fasta
