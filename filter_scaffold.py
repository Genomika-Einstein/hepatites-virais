# SINTAXE: python3 filter_scaffold.py in_path out_path
# Filtra os contigs do FASTA em in_path e guarda o resultado em out_path

import re
from sys import argv

class node():
	def __init__(self, num = -1, cov = -1, length = -1):
		self.fasta = ""
		self.number = num
		self.coverage = cov
		self.length = length

def find_scaffold(path):
# Retorna scaffold (i.e. lista de contigs) presente em um arquivo (scaffold.fasta, do spades)
	file = open(path, 'r')
	lines = file.readlines()
	file.close()
	scaffold = []
	for line in lines:
		regex = r">NODE_([0-9]+)_length_([0-9]+)_cov_([0-9.]+)"
		match = re.search(regex, line)
		if match:
			num = int(match.group(1))
			length = int(match.group(2))
			cov = float(match.group(3))
			scaffold.append(node(num, cov, length))
		scaffold[-1].fasta += line
	return(scaffold)

def filter_scaffold(scaffold):
# Recebe scaffold (i.e. lista de contigs) cru e retorna scaffold filtrado
	length_ratio_cutoff = 0.30 # PARAMETRO
	coverage_cutoff = 10.0 # PARAMETRO
	return [contig for contig in scaffold if contig.coverage > coverage_cutoff and contig.length > length_ratio_cutoff * scaffold[0].length]

def write_scaffold(scaffold, path):
# Recebe scaffold (i.e. lista de contigs) e escreve o fasta correspondente em um arquivo (path)
	contigs = [contig.fasta for contig in scaffold]
	fasta = "".join(contigs)
	file = open(path, 'w')
	file.write(fasta)
	file.close()

def main():
	in_path, out_path = argv[1], argv[2]
	original_scaffold = find_scaffold(in_path)
	new_scaffold = filter_scaffold(original_scaffold)
	write_scaffold(new_scaffold, out_path)

main()
