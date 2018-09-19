import pathlib
import os.path
import sys
import socket

source_path = '/common/home/yusipov_i/Work/mlmg/source/python/Methylation'
source_name = 'infrastructure.hierarchy'

os.system('sbatch ./submit.sh' + ' ' + source_path + ' ' + source_name)
