from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *
from annotation.bop import *

method = Method.manova
host_name = socket.gethostname()
fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
class_type = ClassType.class_a
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big

config = Config(fs_type, db_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type)

dict_bop_cpgs = get_dict_bop_cpgs(config)

dict_cpg_gene = get_dict_cpg_gene(config)
