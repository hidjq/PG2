
# docker build -t idash_t1_pp1_pg1:v1.1 -f Dockerfile  ./
# docker run -itd idash_t1_pp1_pg1:v1.1

FROM ubuntu:20.04
LABEL maintainer="GUO Yunzhe <yunzhe.guo@outlook.com>"

RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list\
 && sed -i s@/security.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list\
 && apt clean\
 && apt-get -y update \
# RUN apt-get install -y software-properties-common && add-apt-repository ppa:deadsnakes/ppa
 && apt-get install -y python3.8\
 && apt-get install -y wget \
 && apt-get install -y unzip \
 && apt-get install -y screen \ 
 && apt-get install -y emacs \ 
 && apt-get install -y libmpfr-dev libmpc-dev \ 
 && rm -f /usr/bin/python3 && ln -s /usr/bin/python3.8 /usr/bin/python3\
 && apt-get install -y python3-pip\
# 清除缓存
 && rm -rf /var/lib/apt/lists/* \
# 修改系统时区
 && rm -rf /etc/localtime && ln -s /usr/share/zoneinfo/Asia/Shanghai /etc/localtime \
# 修改 python 路径
 && rm -f /usr/bin/python \
 && ln -s /usr/bin/python3.8 /usr/bin/python

WORKDIR /root/
COPY . /root/
# 安装python依赖，在这之前需要安装好 gmpy2 所需依赖
RUN pip3 install numpy==1.21.1 -i https://pypi.tuna.tsinghua.edu.cn/simple \
&& pip3 install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple \
&& rm -f requirements.txt && rm Dockerfile && rm -rf __pycache__ 



