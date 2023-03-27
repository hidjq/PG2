
import random 
import pickle
from numba import jit, int32
from tqdm import tqdm
from phe import paillier
from config import random_range_min, random_range_max
from config import Missing, base_list, each_cnt
from config import res_mat_path
from utils.gene import seq_encoding

def generate_key(key_size = 3072):
    public_key, private_key = paillier.generate_paillier_keypair(private_keyring=None,n_length=key_size)
    return public_key, private_key

# 旧方法，不使用
def encrypt_seq(encoded_seq, public_key):
    res = [public_key.encrypt(x) for x in encoded_seq]
    return res 

# 计算两个基因seq的距离，该操作在B方进行
def compare_sequece(enc_seq, origion_seq, cipher_dict,\
    persistence = False, persistence_path = None, pbar = None):
    plain_seq = seq_encoding(origion_seq)
    length = len(plain_seq)
    res = [] 
    for i in range(length):
        enc_index = enc_seq[i]
        enc_base = cipher_dict[enc_index]
        # print(enc_base)
        # print(plain_seq[i])
        if plain_seq[i] == Missing:
            res.append(enc_base * 0)
            # res.append(zero)
        else:
            # 随机化这个 option 比较费时间
            # r = random.uniform(random_range_min, random_range_max)
            # phe.paillier.EncryptedNumber
            compare_value = (enc_base - plain_seq[i]) # * r 
            # print(type(compare_value))
            res.append(compare_value)
    # 打乱顺序防止反推
    random.shuffle(res)
    if pbar is not None:
        pbar.update(1)
    if persistence:
        # write to disk 
        with open(persistence_path, "wb") as f:
            pickle.dump(res, f)
        return 0
    else:
        return res 

# 解密seq计算结果，计算距离分数，不仅是解密操作，该操作在A方运行
def calculate_distance(compare_seq, secret_key):
    length = len(compare_seq)
    # print(length)
    plain = [secret_key.decrypt(x) for x in compare_seq]
    # print("plain")
    # print(plain)
    diff_cnt = 0
    for i in plain:
        if i != 0 and abs(i) < (14-1)*10:
            diff_cnt += 1
    return diff_cnt / length , plain

def calculate_distance_persistence(i,j,secret_key):
    file_path = res_mat_path + "%d_%d.pickle"%(i,j)
    f = open(file_path, 'rb')
    compare_seq = pickle.load(f)
    f.close()

    # 其余计算保持不变
    length = len(compare_seq)
    # print(length)
    plain = [secret_key.decrypt(x) for x in compare_seq]
    print("calculate_distance_persistence.plain")
    print(plain)
    diff_cnt = 0
    for i in plain:
        if i != 0 and abs(i) < (14-1)*10:
            diff_cnt += 1
    return diff_cnt / length , plain


@jit
def calc_cipher_index(base, cipher_cnt ):
    return base * 1000 + cipher_cnt

def  generate_encrypt_hashmap(public_key):    
    # // 13 个 base 含 N
    # base_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1013]
    # // println!("base list: {:?}", base_list);
    # // 建立 hash map 
    hash_map = {}
    pbar = tqdm(total=len(base_list) * each_cnt, desc="Hashmap Encryption", ncols=80)
    for base in base_list:
        enc_base_list = []
        for _ in range(each_cnt):
            enc_base = public_key.encrypt(base)
            enc_base_list.append(enc_base)
            pbar.update(1)
        hash_map[base] = enc_base_list
    pbar.close()
    return hash_map

def ciphertext_randomize(encrypt_hashmap):
    # // 返回 obs_dict 以及 cipher_dict
    # println!("Generating obs and cipher dict...");
    # // 还是建立同样的 base_list 
    # base_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1013]
    # // 开始循环，取数，放在新的 hash map 中
    # rng = rand::thread_rng(); // 用于生成随机数的 rng
    # // 生成需要返回的两个变量，都是 mut 格式
    obs_dict = {}
    cipher_dict = {} 
    # // 开始循环
    r_list = []
    while(len(r_list) < len(base_list) * each_cnt + 10):
        x = random.randint(256, 65536)
        if x not in r_list:
            r_list.append(x)
    for base in base_list:
        for i in range(each_cnt):
            cipher_index = calc_cipher_index(base, i)
            obs_index = r_list.pop()
            cipher_txt = encrypt_hashmap[base][i]
            # // utils::print_type_of(cipher_txt);
            obs_dict[cipher_index] = obs_index
            cipher_dict[obs_index] = cipher_txt
    
    return obs_dict, cipher_dict

# 这个操作用 python 较慢，后期用 rust 代替
# 争取从 15min -> 3min
def encrypt_seqence_replace(obs_dict, v):
    # // 预先生成随机化 rng
    # let mut rng = rand::thread_rng();
    # each_cnt = len(cipher_dict[1])
    res_mat = []
    # print(len(v))
    # /开始迭代
    pbar = tqdm(total=len(v) ,desc="Ciphertxt  Replace", ncols = 80)
    for str_seq in v:
        encode_seq = seq_encoding(str_seq)
        res_seq = []
        for base in encode_seq:
            # // r =  应该在 [0 , cnt) 之间，不含 cnt
            r = random.randint(0, each_cnt-1)
            cipher_index = calc_cipher_index(base , r)
            obs_index = obs_dict[cipher_index]
            # // println!("{}", obs_index);
            res_seq.append(obs_index)
        res_mat.append(res_seq)
        pbar.update(1)
        
    #  // 迭代结束
    pbar.close()
    return res_mat


        
    
