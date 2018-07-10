from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *

num_folds = 10

host_name = socket.gethostname()
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

attributes = get_attributes(config)

cpgs, vals = get_cpg_data(config)

param_names, param_values = get_enet_params(attributes, vals, num_folds)

fn = 'enet_params_cpg.txt'
fn = get_param_path(fs_type, db_type, fn)
save_params(fn, param_names, param_values)

