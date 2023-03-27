import numpy as np 
from Bio import SeqIO
from config import replace_dict, Missing

# 用于算 plain_matrix
def compute_pdistanc_plain(seq1, seq2): 
    score=0
    snp=0
    miss=["N","n","-"]
    for n in range(len(seq1)):
        if seq1[n]!=seq2[n]:
            if seq1[n] in miss or seq2[n] in ['N', 'n', '-']:
                pass
            else:
                snp+=1
    score=snp/len(seq1)
    return format(score,".8f")

# 用于算 plain_matrix
def align_star_multiple(seqs):
    scores = np.zeros((len(seqs), len(seqs)))
    for x in range(len(seqs)):
        for y in range(len(seqs)):
            if x!=y:
                scores[x,y]=compute_pdistanc_plain(seqs[x], seqs[y])
    return scores

def read_seq_from_disk(fasta_path = "../data/L1/L1-aln.fasta",\
    percent = 1, backward = False):
    Seq=[]
    for seq_record in SeqIO.parse(fasta_path,"fasta"): 
        Seq.append(str(seq_record.seq))
    length = len(Seq)
    ret_length = int(length * percent)
    # all gene sequences are stored in seq 
    if backward is False:
        return Seq[:ret_length]
    else:
        return Seq[(-1)*ret_length:]

# ID 需要作为 nj 数算法的参数
def read_id_from_disk(fasta_path = "../data/L1/L1-aln.fasta",\
    percent = 1, backward = False):
    ID=[]
    for seq_record in SeqIO.parse(fasta_path, "fasta"): 
        ID.append(seq_record.id)
    length = len(ID)
    ret_length = int(length * percent)
    # all gene sequences are stored in seq 
    if backward is False:
        return ID[:ret_length]
    else:
        return ID[(-1)*ret_length:]

# 针对单条基因 seq 
def seq_encoding(input_seq):
    res = []
    for e in input_seq:
        res.append(replace_dict[e])
    return res 

def compute_pdistance(seq1, seq2): 
    score=0
    snp=0
    miss=["N","n","-", Missing]
    for n in range(len(seq1)):
        if seq1[n]!=seq2[n]:
            if seq1[n] in miss or seq2[n] in miss:
                pass
            else:
                snp+=1
    score=snp/len(seq1)
    print(snp)
    return format(score,".8f")