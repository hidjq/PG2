
import time 
import os

# 删除上次运行的文件
def clean_up():
    cmd_pickle = '''find ./ -name '*.pickle' -type f -exec rm {} \;  '''
    cmd_json   = '''find ./ -name '*.json' -type f  -exec rm {} \;  '''
    # cmd_tree   = '''find ./ -name '*.tree' -type f  -exec rm {} \;  '''

    os.system(cmd_pickle)
    os.system(cmd_json)
    # os.system(cmd_tree)

def create_dir():
    os.system('mkdir -p B_to_A/res_mat/')
    os.system('mkdir -p A_to_B/')
    os.system('mkdir -p A_local/')
    os.system('mkdir -p output/')
    pass


def note_process(process,party = 'A',  Done = False):
    assert party in ['A', 'B']
    if party == 'A': # in yellow 
        party_str = '\033[0;33;40mParty A\033[0m'
    else: # 
        party_str = '\033[0;36;40mParty B\033[0m'
    done_str = '\033[0;32;40mDone.\033[0m\n' if Done else "Start."
    print(time.ctime(), "- %s: " % party_str ,process, done_str)

def get_dir_size(dir):
    size = 0
    for root, dirs, files in os.walk(dir):
        size += sum([os.path.getsize(os.path.join(root, name)) for name in files])
    return size / 1024 / 1024 # 返回 MB
