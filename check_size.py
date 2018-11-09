from sys import argv
file=open(argv[1])
for line in file.readlines():
	if line[0]!='>':
		L = line.strip()
		print(len(L), end='\t')
		print(100*len([c for c in L if c!='N'])/len(L))
	else:
		print(line.strip(), end='\t')
file.close()
