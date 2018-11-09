for i in $(cat samples.txt)
do
	clear
	clear
	if [ -e ${i}/${i}_scaffold.fasta ]
	then
		python3 fasta_reader.py ${i}/${i}_scaffold.fasta
		cat ${i}/${i}_scaffold.fasta
	fi
	if [ -e ${i}/${i}_blast_filtered.txt ]
	then
		cat ${i}/${i}_blast_filtered.txt
		read var
	fi
done
