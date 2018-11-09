sample="$1" # sample name (string)
data="/home/ubuntu/bioinfo/data"
echo "Iniciando busca..."
blastn -query $data/contigs/${sample}_contigs.fasta -remote -db refseq_genomic_human -outfmt 6 -out $data/hgblast/${sample}_hgblast.csv \
-word_size 28 -evalue 1 -perc_identity 98 -qcov_hsp_perc 95
echo "Busca terminada."
