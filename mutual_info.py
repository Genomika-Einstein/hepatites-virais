from sys import argv
from sklearn.metrics import mutual_info_score
from numpy import log
from matplotlib import pyplot as plt

def letter2int(s):
	t = []
	for base in s:
		if base.lower() == 'a':
			t.append(1)
		elif base.lower() == 't':
			t.append(2)
		elif base.lower() == 'c':
			t.append(3)
		elif base.lower() == 'g':
			t.append(4)
		else:
			t.append(0)
	return t

def get_fasta(path):
	file = open(path, 'r')
	lines = [line.strip() for line in file.readlines()]
	fasta = {}
	for line in lines:
		if line[0] == '>':
			header = line
			fasta[header] = ""
		else:
			fasta[header] += line
	return fasta

def genotype(header):
	return ord(header.split('_')[-1][0].lower())-ord("a")

def distribution(fasta):
	# P[i] é a matrix p tal que p[x][y] é a Prob(genotipo==x and seq[i]==y)
	P = [[[0 0 0 0 0] for i in range(10)] for j in len(fasta[0])]
	for header in fasta:
		seq = letter2int(fasta[header])
		for pos in range(len(seq)):
			P[pos][genotype(header)][seq[pos]] +=1
	return P

def main():
	fasta = get_fasta(argv[1])
	L = len(fasta[list(fasta.keys())[0])
	P = distribution(fasta)
	MI = [mutual_info_score(P[i])/log(2) for i in range(L)]
	print(MI)
	plt.plot(range(1,len(MI)+1), MI)
	savefig("mutual_info.png")

main()