from sys import argv

CODE = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
	}

def translate(DNA):
	peptide = ""
	num_amino = len(DNA)//3 # numero de aminoacidos
	for pos in range(0, 3*num_amino, 3):
		peptide += CODE[DNA[pos:pos+3]]
	if peptide == '':
		return '.'
	else:
		return peptide

def find(pos, track):
# dada uma posição (1-based) no genoma em coordenadas absolutas, 
# encontra a posição (0-based) em relação ao começo do gene 
# (dado como uma lista de tuplas (ini, fim)).
	pos = pos - 1 # 1- to 0-based
	pos_rel = 0
	for (ini, fim) in track:
		if ini <= pos and pos < fim:
			pos_rel += pos - ini
			return pos_rel
		else:
			pos_rel += fim - ini
	raise ValueError("position "+str(pos+1)+"not found in track "+track)

def isin(pos, track):
	pos = pos - 1 # 1- to 0-based
	for (ini, fim) in track:
		if ini <= pos and pos < fim:
			return True
	return False

def get_vcf(path):
	file = open(path, 'r')
	vcf = {} # vcf = {chrom: [(pos, ref, alt)]}
	for line in file.readlines():
		if line[0] != '#':
			L = line.strip().split('\t')
			chrom, pos, ref, alt = L[0], int(L[1]), L[3], L[4]
			if chrom not in vcf:
				vcf[chrom] = []
			vcf[chrom].append( (pos, ref, alt) )
	file.close()
	return vcf

def get_fasta(path): # fasta = {track: seq}
	file = open(path, 'r')
	fasta = {}
	for line in file.readlines():
		if line[0] == '>':
			header = line[1:-1]
			fasta[header] = ""
		else:
			fasta[header] += line.strip()
	file.close()
	return fasta

def get_bed(path): # bed = {chrom: {track: [(ini, fim)]}}
	file = open(path, 'r')
	bed = {}
	for line in file.readlines():
		L = line.strip().split('\t')
		chrom, ini, fim, track = L[0], int(L[1]), int(L[2]), L[3]
		if chrom not in bed:
			bed[chrom] = {}
		if track not in bed[chrom]:
			bed[chrom][track] = []
		pair = (ini, fim)
		bed[chrom][track].append(pair)
	file.close()
	return bed

def find_mutation(seq1, seq2):
	# dadas duas sequencias seq1 != seq2, retorna (i,j)
	# tais que a a diferença entre seq1 e seq2 esteja toda
	# localizada em seq1[i:j+1].
	# OBS: i>0, j<0
	i = 0
	while seq1[i] == seq2[i]:
		i = i+1
	j = -1
	while seq1[j] == seq2[j]:
		j = j-1
	return (i,j)

def output(text, path):
	file = open(path, 'a')
	file.write(text+'\n')
	file.close()

# variaveis p/ contagem
snp = 0
mnp = 0
frameshift = 0
syn = 0

# abrindo arquivos
vcf = get_vcf(argv[1])
fasta = get_fasta(argv[2])
bed = get_bed(argv[3])
for path in ['.bcp.txt', '.frm.txt', '.ptm.txt', '.syn.txt']:
	file = open(argv[4]+path, 'w')
	file.write('')
	file.close()

for chrom in bed:
	for line in vcf[chrom]:
		pos, ref, alt = line[0], line[1], line[2]
		for track in bed[chrom]:
			if track == 'Pol' and isin(pos, bed[chrom]['RT']):
				continue # avoid double reporting RT mutations
			if isin(pos, bed[chrom][track]):
				pos_rel = find(pos, bed[chrom][track])
				dna_wt = fasta[track]
				dna_mut = dna_wt[:pos_rel] + alt + dna_wt[pos_rel+len(ref):]
				prot_wt = translate(dna_wt)
				prot_mut = translate(dna_mut)
				if track == 'BCP': # promotores e enhancers, nt EcoRI
					output(ref + str(pos) + alt, argv[4]+'.bcp.txt')
				elif (len(dna_wt)-len(dna_mut))%3 != 0: # Frameshift, nt ORF
					output(track.lower()+ref+str(pos_rel)+alt, argv[4]+'.frm.txt')
					frameshift += 1
				elif prot_wt != prot_mut: # mutação s/ frameshift, aa ORF
					i, j = find_mutation(prot_wt, prot_mut)
					output(track.lower()+prot_wt[i:j+1] + str(pos_rel//3+1) + prot_mut[i:j+1], argv[4]+'.ptm.txt')
					if len(prot_mut[i:j+1]) == 1:
						snp += 1
					else:
						mnp += 1
				else: # prot_wt == prot_mut, mutação sinônima
					output(track.lower()+ref+str(pos_rel+1)+alt, argv[4]+'.syn.txt')
					syn += 1
print('single\t',snp,'\ncomplex\t',mnp,'\nsynonymous\t',syn,'\nframeshift\t',frameshift)
