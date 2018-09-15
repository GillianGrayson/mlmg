import pathlib
import os.path


source_path = '/common/home/yusipov_i/Work/mlmg/source/python/Methylation/infrastructure'
source_name = 'hierarchy.py'

os.system('sbatch submit.sh' + ' ' + source_path + ' ' + source_name)
