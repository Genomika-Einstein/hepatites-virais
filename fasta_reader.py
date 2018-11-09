from sys import argv

def get_fasta(path):
	file = open(path)
	hash = {}
	for line in file.readlines():
		if line[0] == ">":
			header = line[1:-1]
			hash[header] = ""
		else:
			hash[header] += line[:-1]
	return hash

fasta = get_fasta(argv[1])
for seq in fasta:
	print(seq, len(fasta[seq]), len([base for base in fasta[seq] if base!='N']))
