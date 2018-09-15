import pathlib
import os.path
import sys

sys.path.append('/common/home/yusipov_i/Work/mlmg/source/python/Methylation')
for p in sys.path:
    print(p)
source_path = '/common/home/yusipov_i/Work/mlmg/source/python/Methylation/infrastructure'
source_name = 'hierarchy.py'

os.system('sbatch ./submit.sh' + ' ' + source_path + ' ' + source_name)
