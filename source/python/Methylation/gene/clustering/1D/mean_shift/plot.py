import scipy.stats as stats
from config import *
from infrastructure.load import *
from infrastructure.save import *
from method.method import *
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.cluster import MeanShift, estimate_bandwidth

num_top_genes = 100
gene_id = 0

method = Method.clustering
clustering_type = ClusteringType.mean_shift
gd_type = GeneDataType.mean
geo_type = GeoType.islands_shores
host_name = socket.gethostname()
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
config = Config(fs_type, db_type, geo_type=geo_type, clustering_type=clustering_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type=geo_type, clustering_type=clustering_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type=geo_type, clustering_type=clustering_type)

attributes = get_attributes(config)

gene_names = get_top_gene_names(config, gd_type, method, num_top_genes)
gene_vals = get_top_gene_vals(config, gd_type, gene_names)

name = gene_names[gene_id]
vals = np.array(gene_vals[gene_id]).reshape(-1, 1)
bandwidth = estimate_bandwidth(vals)
labels = MeanShift(bandwidth=bandwidth).fit_predict(vals)

plt.scatter(list(vals), labels, c=labels, s=100)
plt.title(name)
plt.rcParams.update({'font.size': 30})
plt.show()
