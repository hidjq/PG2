from mpyc.runtime import mpc

import pickle
from tqdm import tqdm
import numpy as np 
import json 
from utils.nj_tree import DistanceMatrix
from utils.nj_tree import nj

from utils.mycrypto import calculate_distance_persistence,compare_sequece, encrypt_seqence_replace
from utils.gene import read_id_from_disk, read_seq_from_disk, compute_pdistance, seq_encoding, align_star_multiple
from config import res_mat_path, pool_size

# A ,B运行
def encryption(fasta_path = "../data/L1/L1-aln.fasta", percent = .5, name = 'a'):
    v = read_seq_from_disk(fasta_path, percent = percent)
    res_mat = encrypt_seqence_replace(v) 
    # print(time.ctime())
    # # 开始持久化
    if(name == 'A'):
        with open("A_to_B/enc_mat_A.pickle", 'wb') as f :
            pickle.dump(res_mat, f)
    else:
        with open("B_to_A/enc_mat_B.pickle", 'wb') as f :
            pickle.dump(res_mat, f)
    # 返回，也可以不返回
    return res_mat

# A, B 都需要用到
def calculation_plain(fasta_path , percent, backward=False):
    plain_ID  = read_id_from_disk(fasta_path, percent=percent, backward = backward)
    plain_seq = read_seq_from_disk(fasta_path, percent=percent, backward = backward) 
    score_matrix=align_star_multiple(plain_seq) 
    return score_matrix, plain_ID

# 实际距离计算，B方执行
def calculation_enc(fasta_path_B = "../data/L1/L1-aln.fasta", percent_B = .5):
    # 先反序列化
    f_mat = open("A_to_B/enc_mat_A.pickle", 'rb')
    enc_mat_A = pickle.load(f_mat)
    f_mat.close()  
    
    # 之后要改为 True
    v = read_seq_from_disk(fasta_path_B, percent=percent_B, backward=True) 
    enc_mat_B = encrypt_seqence_replace(v)
    
    # length_A = 5  # len(enc_mat)
    # length_B = 5  # len(plain_mat)

    length_A = len(enc_mat_A)
    length_B = len(enc_mat_B)
    # till now the enc_matA and enc_matB a two part
    # print(length_A, length_B)
    
    # 改成用角标取放元素，便于多线程处理
    global res_mat
    res_mat = np.zeros(length_A * length_B).reshape(length_A, length_B).tolist()

    # 增加多线程，貌似并没有变快？
    # pool = Pool(pool_size)
    # results = []
    pbar = tqdm(total=length_A * length_B, desc="Distance Calculation(HE)", ncols = 80 )
    for i in range(length_A):
        for j in range(length_B):
            pickle_path = res_mat_path + "%d_%d.pickle"%(i,j)
            # pool.apply_async(compare_sequece, (enc_mat[i],plain_mat[j],cipher_dict,True,  pickle_path, pbar))
            # results.append(result)
            # 以下为单线程
            seq_res = compare_sequece(enc_mat_A[i],enc_mat_B[j], persistence=True, persistence_path= pickle_path)
            pbar.update(1)
    # [result.wait() for result in results]
    # pool.close()
    # pool.join()

    pbar.close()
    return length_A, length_B
    
# A 方执行
def decryption(length_A, length_B):
    # 反序列化得到私钥
    # global matrix
    matrix =[[[]for i in range(length_A)]for i in range(length_B)]

    # 这里期望可以多线程
    pbar = tqdm(total=length_A * length_B,desc="Distance  Decryption(HE)", ncols=80)
    for i in range(length_A):
        for j in range(length_B):
            pbar.update(1)
            dist, _ = calculate_distance_persistence(i,j)
            # print(mpc.run(mpc.output(dist)))
            matrix[i][j] = dist #dist is cipher
    pbar.close()
    # matrix 为所求
    for i in matrix:
        print(mpc.run(mpc.output(i)))
    return matrix

# A 方执行建树
def calculate_NJ_tree(score_matrix, ID, output_path):
    dm = DistanceMatrix(score_matrix, ID)
    tree=nj(dm) 
    # 输出为所求
    tree.write(output_path)
    
def compose_score_matrix(matrix_enc, matrix_A, matrix_B):
    length_A = len(matrix_A)
    length_B = len(matrix_B)
    length = length_A + length_B
    #  创建矩阵
    score_matrix = np.zeros(length * length).reshape(length, length)
    # 首先把两个对角放进去
    for i in range(length_A):
        for j in range(length_A):
            score_matrix[i][j] = matrix_A[i][j]
    for i in range(length_B):
        for j in range(length_B):
            score_matrix[length_A+i][length_A+j] = matrix_B[i][j]
    # 然后放同态计算得到的矩阵数值
    for i in range(length_A):
        for j in range(length_B):
            score_matrix[i][length_A+j] = matrix_enc[i][j]
            score_matrix[length_A+j][i] = matrix_enc[i][j]

    return score_matrix
    