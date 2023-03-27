
import pickle
from tqdm import tqdm
import numpy as np 
import json 
from utils.nj_tree import DistanceMatrix
from utils.nj_tree import nj

from utils.crypto import calculate_distance_persistence, generate_key, generate_encrypt_hashmap, ciphertext_randomize, encrypt_seqence_replace
from utils.crypto import compare_sequece, calculate_distance
from utils.gene import read_id_from_disk, read_seq_from_disk, compute_pdistance, seq_encoding, align_star_multiple
from config import res_mat_path, pool_size

# A 方运行
def key_gen(key_size):
    public_key, secret_key = generate_key(key_size=key_size)
    # dump 到硬盘
    # 这两个路径还可以修改
    # public key 
    f_pk = open("A_local/pk.pickle", 'wb')
    pickle.dump(public_key, f_pk)
    f_pk.close()
    # secret key 
    f_sk = open("A_local/sk.pickle", 'wb')
    pickle.dump(secret_key, f_sk)
    f_sk.close()
    return public_key

# A 方运行
def encryption(public_key,fasta_path = "../data/L1/L1-aln.fasta", percent = .5):
    # 首先由于仍然在A方运行，public key 作为参数传入
    # print(time.ctime())
    encrypt_hashmap = generate_encrypt_hashmap(public_key)
    # print(time.ctime())
    obs_dict, cipher_dict = ciphertext_randomize(encrypt_hashmap)
    # obs_dict 使用 json 持久化
    # cipher_dict 可使用其他格式
    # print(time.ctime())
    # 一下内容可以由 rust 代替 
    # print("reading fasta")
    v = read_seq_from_disk(fasta_path, percent = percent)
    # print(time.ctime())
    res_mat = encrypt_seqence_replace(obs_dict, v) 
    # res_mat 使用 json 持久化，方便与 rust 混编
    # print(time.ctime())

    # 开始持久化
    with open("A_to_B/obs_dict.json", 'w') as f:
        json.dump(obs_dict, f)
    with open("A_to_B/cipher_dict.pickle", 'wb') as f:
        pickle.dump(cipher_dict, f)
    with open("A_to_B/enc_mat.json", 'w') as f :
        json.dump(res_mat, f)
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
    f = open("A_to_B/cipher_dict.pickle", 'rb')
    cipher_dict = pickle.load(f)
    f.close()

    f_mat = open("A_to_B/enc_mat.json", 'rb')
    enc_mat = json.load(f_mat)
    f_mat.close()  
    
    # 之后要改为 True
    plain_mat = read_seq_from_disk(fasta_path_B, percent=percent_B, backward=True) 
    
    # length_A = 5  # len(enc_mat)
    # length_B = 5  # len(plain_mat)

    length_A = len(enc_mat)
    length_B = len(plain_mat)

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
            seq_res = compare_sequece(enc_mat[i],plain_mat[j],cipher_dict, persistence=True, persistence_path= pickle_path)
            pbar.update(1)
    # [result.wait() for result in results]
    # pool.close()
    # pool.join()

    pbar.close()
    return length_A, length_B
    
# A 方执行
def decryption(length_A, length_B):
    # 反序列化得到私钥
    f_sk = open("A_local/sk.pickle", 'rb')
    sk = pickle.load(f_sk)
    f_sk.close()

    # global matrix
    matrix = np.zeros(length_A * length_B).reshape(length_A, length_B)

    # 这里期望可以多线程
    pbar = tqdm(total=length_A * length_B,desc="Distance  Decryption(HE)", ncols=80)
    for i in range(length_A):
        for j in range(length_B):
            pbar.update(1)
            # dist, _ = calculate_distance(res_mat[i][j], sk)
            dist, _ = calculate_distance_persistence(i,j, sk)
            matrix[i][j] = dist
    pbar.close()
    # matrix 为所求
    # print(matrix)
    return matrix.tolist()

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
    