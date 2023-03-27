

import time 
import pickle
import sys
import os 
from utils.proc import calculate_NJ_tree, calculation_plain, compose_score_matrix, key_gen, encryption, calculation_enc, decryption
from utils.utils import clean_up, note_process, create_dir, get_dir_size
from utils.extract_data import extraction_proc
from mpyc.runtime import mpc

from config import percent_A, percent_B
assert(percent_A + percent_B == 1)

# 按需求解压数据，如果以前解压过，就不再次解压
extraction_proc()

# 处理命令行参数
try:
    index_str = sys.argv[1]
    if index_str in ['L1', '1']:
        index = 1
    elif index_str in ['L2', '2']:
        index = 2
    elif index_str in ['L3', '3']:
        index = 3
    else:
        index = 0
except IndexError:
    index = 0

fasta_path = "./data/L%d/L%d-aln.fasta" % (index,index)
print("==" * 16,"\nUse Data: %s" %fasta_path)
print("==" * 16, "\n")


# 预备工作
clean_up()
create_dir()

# A
st_time = time.time()
note_process("Key Generation (3072-bit)")
public_key = key_gen(key_size = 3072)
note_process("Key Generation (3072-bit)", Done=True)
ed_time = time.time()
print("Time of Key Generation: %g seconds\n" % (ed_time - st_time))

st_time = time.time()
note_process("Gene Data Encryption")
enc_mat = encryption(public_key, fasta_path=fasta_path, percent=percent_A)
note_process("Gene Data Encryption", Done=True)
ed_time = time.time()
print("Time of Gene Data Encryption: %g seconds\n" % (ed_time - st_time))


# B
st_time = time.time()
note_process("Distance Calculation(in Secret)", party='B')
length_A, length_B = calculation_enc(fasta_path_B = fasta_path, percent_B=percent_B)
note_process("Distance Calculation(in Secret)", party='B', Done = True)

note_process("Score_Matrix(Plain) Calculation", party='B')
score_matrix_B, ID_B = calculation_plain(fasta_path,percent=percent_B, backward=True)
note_process("Score_Matrix(Plain) Calculation", party='B', Done=True)
ed_time = time.time()
print("Time of Ciphertext Calculation: %g seconds\n" % (ed_time - st_time))

# A
st_time = time.time()
note_process("Result Decryption")
score_matrix_enc = decryption(length_A, length_B)
note_process("Result Decryption", Done=True)
ed_time = time.time()
print("Time of Decryption: %g seconds\n" % (ed_time - st_time))

print("Communication Size A -> B: %20.4f MB." % get_dir_size("A_to_B/"))
print("Communication Size B -> A: %20.4f MB.\n" % get_dir_size("B_to_A/"))


# 计算tree (A方)
note_process("Prepare dm matrix")
score_matrix_A, ID_A = calculation_plain(fasta_path,percent=percent_A, backward=False)
final_score_matrix = compose_score_matrix(score_matrix_enc, score_matrix_A, score_matrix_B)
note_process("Prepare dm matrix", Done=True)

note_process("Build NJ Tree")
calculate_NJ_tree(final_score_matrix, ID_A + ID_B, output_path = "output/nj_tree_L%d.tree"%index)
note_process("Build NJ Tree", Done= True)
print("NJ Tree path:", "output/nj_tree_L%d.tree"%index)

# 输出 score matrix
f = open("output/final_score_matrix_L%d.pickle"%index, 'wb')
pickle.dump(final_score_matrix,f)
f.close()

# 输出 病毒序列ID到 txt文件
# score matrix 行列顺序相同，对角线距离数值为0

ID_list  =ID_A + ID_B
f = open('output/id_list_row_L%d.txt' % index, 'w')
f.writelines("\n".join(ID_list))
f.close()

f = open('output/id_list_col_L%d.txt' % index, 'w')
f.writelines("\n".join(ID_list))
f.close()

os.system("rm -rf __pycache__") # 删除影响观感的目录

