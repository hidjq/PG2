# coding=utf-8

import zipfile
import os 
import sys
import os

if sys.version_info < (3, 0):
    print("Please use Python 3!")
    print("Exit.")
    sys.exit(0)


def extraction_proc():

    if os.path.exists("data/L0/L0-aln.fasta"):
        return 
    

    UNZIP_PASSWD = '2021JYSD$wppcc'
    ZIP_PATH = "202111-data.zip"

    print("Unzip Files Start.")
    # 解压文件
    # 解压后文件会通过 .gitignore 忽略，提交git时候不会提交至仓库，避免仓库体积过大
    zip_file = zipfile.ZipFile(ZIP_PATH)
    zip_list = zip_file.namelist()
    for f in zip_list:
        zip_file.extract(f, './data/',pwd=UNZIP_PASSWD.encode("utf-8"))
    zip_file.close()
    

    # print("Moving Files start.")
    # 移动文件至 data/ 并删除无关目录
    os.system("mv data/202111*/* data") # 移动L1/L2/L3目录到 ./data/
    os.system("find data/ -type d -empty -delete") # 删除空目录
    # 建立一个更短的 L0 数据集
    # print("Moving Files Done.")

    # print("Creating L0 Data start.")
    os.system("mkdir -p data/L0")
    os.system("head -n 20 data/L1/L1-aln.fasta > data/L0/L0-aln.fasta")
    # print("Creating L0 Data Done.")

    # os.system("mv 202111-data.zip .202111-data.zip")
    print("Unzip Files \033[0;32;40mDone.\033[0m\n")
    
