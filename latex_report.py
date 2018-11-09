from sys import argv

def HGcontamination(path):
	tex = "\\section{Contaminação por DNA humano}\n"
	file = open(path, 'r')
	L = file.readlines()
	if len(L) == 0:
		tex += "Não há evidências de contaminação da amostra.\n"
	else:
		tex += "Há evidência de contaminação humana. Segue relatório do BLAST:\\\\\n"
		tex += "\\begin{tabular}"
		tex += "Contig & Referência & Identidade (\\%) & Tamanho (nt) & Gaps & Início (contig) & Fim (contig) & Início (ref.) & Fim (ref.) & E-value & Bit score"
		tex += "\\csvreader"
		tex += "\\end{tabular}"


tex =  ""
tex += "\\usepackage{csvsimple}"
