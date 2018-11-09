# igv_caller.sh
# configura o Xvfb e chama o igv.
# SINTAXE: igv_caller.sh -bam input.bam -batch batch_script.txt -g genome.fasta

igv_dir="/home/ubuntu/bioinfo/app/IGV_2.4.14"
ref="/home/ubuntu/bioinfo/reference/HBV_ref/HBV_ref.fasta"
bed="/home/ubuntu/bioinfo/reference/HBV_ref/HBV_ref.bed"
data="/home/ubuntu/bioinfo/data"
batch="igv_batch.txt"

while [ "$1" != "" ]
        do
        case $1 in
                "-sample")
                        shift
                        sample="$1"
                        shift
                        ;;
                "-ref")
                        shift
                        ref="$1"
                        shift
                        ;;
        esac
        done

# make sure Xvfb is not up
Xvfb_pid=$(pgrep Xvfb)
if [ ! -z "$Xvfb_pid" ]
	then kill -9 $Xvfb_pid
fi
# launch Xvfb
Xvfb :0 -nolisten tcp &
export DISPLAY=:0

# make batch file
echo new > $batch
echo load $data/BAM/$sample.bam >> $batch
echo load $bed >> $batch
echo expand $(basename $bed) >> $batch
echo snapshotDirectory $data/figures/igv >> $batch
echo snapshot ${sample}_igv.png >> $batch
echo exit >> $batch

# start igv
bash $igv_dir/igv.sh $data/BAM/$sample.bam -g $ref -b $batch
