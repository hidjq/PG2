# 赛题相关可配置参数
# 甲乙两方数据占比，须相加为 1
percent_A = .5
percent_B = .5

# 其中 percent_A * 样本总数量 = M 
# 其中 percent_B * 样本总数量 = N 

# 基因编码自定义区域
A = 1
G = 2
C = 3
T = 4
Y = 5 
M = 6
R = 7
V = 8
W = 9
B = 10
K = 11
S = 12
H = 13
D = 14
Missing = 1015
base_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14, Missing]

replace_dict = {
    "A":A, "G":G, "C":C, "T":T,"D": D, # 正常基因
    "Y":Y, "M": M, "R": R, "V": V, "W": W,"B": B, "K": K,"S": S, "H": H,
    "N":Missing, "-":Missing # 缺失基因
}
# Hashmap 每个碱基随机化密文数量配置
each_cnt = 42

# 程序辅助变量
zero_append_probability = .2
pool_size = 12
res_mat_path = "B_to_A/res_mat/"
random_range_min = 1
random_range_max = 10

