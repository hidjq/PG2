from mpyc.runtime import mpc

import random 
import pickle
from numba import jit, int32
from tqdm import tqdm
from phe import paillier
from config import random_range_min, random_range_max
from config import Missing, base_list, each_cnt
from config import res_mat_path
from utils.gene import seq_encoding

secint = mpc.SecInt()
secfxp = mpc.SecFxp()

def encrypt_seqence_replace(v):
    res_mat = []
    # /开始迭代
    pbar = tqdm(total=len(v) ,desc="Ciphertxt  Replace", ncols = 80)
    for str_seq in v:
        encode_seq = seq_encoding(str_seq)
        res_seq = []
        for base in encode_seq:            
            ciper_index = secint(base)
            res_seq.append(ciper_index)
        res_mat.append(res_seq)
        pbar.update(1)
        
    #  // 迭代结束
    pbar.close()
    # for i in res_mat:
    #     for j in i:
    #         print(mpc.run(mpc.output(j)),end='')
    #         print(' ',end='')
    return res_mat

# 计算两个基因seq的距离，该操作在mpc进行
def compare_sequece(enc_seq_A, enc_seq_B,\
    persistence = False, persistence_path = None, pbar = None):
    length = len(enc_seq_B)
    res = [] 
    # print('seq A plain:')
    # print(mpc.run(mpc.output(enc_seq_A)))
    # print('seq B plain:')
    # print(mpc.run(mpc.output(enc_seq_B)))
    for i in range(length):

        # print('calculate scp = ')
        # print(mpc.run(mpc.output(mpc.if_else(mpc.eq(enc_seq_A[i], Missing),0 ,
        # mpc.if_else(mpc.eq(enc_seq_B[i],Missing),0 ,enc_seq_A[i]-enc_seq_B[i])))))
        k = secint(mpc.run(mpc.output(mpc.if_else(mpc.eq(enc_seq_A[i], Missing),0 ,
            mpc.if_else(mpc.eq(enc_seq_B[i],Missing),0 ,enc_seq_A[i]-enc_seq_B[i])))))
        res.append(k)
    # 打乱顺序防止反推
    # print('res = ')
    # print(mpc.run(mpc.output(res)))
    if pbar is not None:
        pbar.update(1)
    if persistence:
        # write to disk 
        with open(persistence_path, "wb") as f:
            pickle.dump(res, f)
        return 0
    else:
        return res 

def calculate_distance_persistence(i,j):
    file_path = res_mat_path + "%d_%d.pickle"%(i,j)
    f = open(file_path, 'rb')
    compare_seq = pickle.load(f)
    f.close()

    print('compare_seq = ')
    print(mpc.run(mpc.output(compare_seq)))
    # 其余计算保持不变
    length = len(compare_seq)
    # print(length)
    # print(plain)
    diff_cnt = secint(0)
    for i in compare_seq:
        diff_cnt = diff_cnt + ~mpc.eq(i,0)
    # print('diff_cnt = ')
    # print(mpc.run(mpc.output(diff_cnt)))

    flt = mpc.convert(diff_cnt, secfxp)
    # print('flt = ')
    # print(mpc.run(mpc.output(flt)))
    # print(mpc.run(mpc.output(mpc.div(flt, length))))

    return mpc.div(flt, length) , compare_seq