# NODE_1_length_3087_cov_187.783108     gi|22135701|gb|AY090455.1|HBV_subgenotipo_F2a   98.39   1678    27      0       1407    3084    1       1678    0.0     2950
# 0       1       2       3       4         5        6       7     8       9     10      11		12
# qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore,	strand

from sys import argv

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0

def max_key(jisho):
    res = -1
    for key in jisho:
        if key > res:
            res = key
    return res

def genotype(hit):
    return hit[1].split('_')[-1]

def matrix_to_blast_txt(matrix, txt_path):
    file = open(txt_path, 'a')
    for row in matrix:
        file.write('\t'.join([str(x) for x in row]) + '\n')
    file.close()

def blast_txt_to_matrix(txt_path):
    file = open(txt_path, 'r')
    txt = file.readlines()
    file.close()
    matrix = []
    for line in txt:
        hit = line.split('\t')
        hit = [hit[0], hit[1], float(hit[2])] + [int(x) for x in hit[3:10]] + [float(hit[10]), float(hit[11])]
        matrix.append(hit)
    return matrix

def fasta_to_hash(fasta_path):
    file = open(fasta_path)
    fasta = file.readlines()
    file.close()
    hash = {}
    for line in fasta:
        if line[0] == ">":
            header = line[1:-1]
            hash[header] = ""
        else:
            hash[header] += line[:-1]
    return hash

def separate_by_ref(blast_table):
    hash = {}
    seen_genotypes = []
    for hit in blast_table:
        if hit[1] in hash:
            hash[hit[1]].append(hit)
        if genotype(hit) not in seen_genotypes:
            seen_genotypes.append(genotype(hit))
            hash[hit[1]] = [hit]
    return hash

def remove_intersec(blast_table):
    n = len(blast_table)
    qini = [hit[6] for hit in blast_table]
    qfim = [hit[7] for hit in blast_table]
    sini = [hit[8] for hit in blast_table]
    sfim = [hit[9] for hit in blast_table]
    sense = [1 if sini[i] < sfim[i] else -1 for i in range(n)]
    for worse in reversed(range(n)):
        for better in range(worse):
            if qini[better] <= qini[worse] and qini[worse] <= qfim[better]:
                sini[worse] += sense[worse]*(qfim[better] + 1 - qini[worse])
                qini[worse] = qfim[better] + 1
            if qini[better] <= qfim[worse] and qfim[worse] <= qfim[better]:
                sfim[worse] += sense[worse]*(qini[better] - 1 - qfim[worse])
                qfim[worse] = qini[better] - 1
    redundantes = [i for i in range(n) if qini[i] >= qfim[i]]
    new_table = []
    for i in range(n):
        if i not in redundantes:
            old_row = blast_table[i]
            new_row = old_row[:6] + [qini[i], qfim[i], sini[i], sfim[i]] + old_row[10:]
            new_table.append(new_row)
    return new_table

def previous_assembly(blast_table, contigs_hash):
    sequence = {}
    for hit in reversed(blast_table):
        contig = contigs_hash[hit[0]]
        for contig_pos in range(hit[6], hit[7]+1):
            subject_pos = contig_pos + hit[8] - hit[6]
            try:
                sequence[subject_pos] = contig[contig_pos]
            except:
                continue
    string = ""
    for i in range(max_key(sequence)):
        if i in sequence:
            string += sequence[i].upper()
        else:
            string += "N"
    return string

def assembly(blast_table, contigs_hash):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    size = {'A':3221, 'B':3215, 'C':3215, 'D':3182, 'E':3212, 'F':3215, 'G':3248, 'H':3215, 'I':3215, 'J':3182}
    sequence = {}
    for hit in reversed(blast_table):
        gtype = genotype(hit)
        contig = contigs_hash[hit[0]]
        contig_pos = hit[6]
        subject_pos = hit[8]
        sense = sign(hit[9]-hit[8])
        while contig_pos <= hit[7]:
            if sense>0:
                sequence[subject_pos-1] = contig[contig_pos-1]
            elif sense<0:
                sequence[subject_pos-1] = complement[contig[contig_pos-1].upper()]
            contig_pos += 1
            subject_pos += sense
    string = ""
    for i in range(max_key(sequence)):
        if i in sequence:
            string += sequence[i].upper()
        else:
            string += "N"
    return string + (size[gtype[0]]-len(string))*"N"

def main():
    contigs_hash = fasta_to_hash(argv[2])
    blast_table = blast_txt_to_matrix(argv[1])
    try:
        title = argv[3]
    except:
        title = ""
    blast_table_hash = separate_by_ref(blast_table)
    # apaga o blast.txt original
    f = open(argv[1], 'w')
    f.write('')
    f.close()
    # pronto
    for reference in blast_table_hash:
        blast_table = blast_table_hash[reference]
        blast_table = remove_intersec(blast_table)
        matrix_to_blast_txt(blast_table, argv[1])
        print(">" + title + "/" + reference)
        print(assembly(blast_table, contigs_hash))

main()
