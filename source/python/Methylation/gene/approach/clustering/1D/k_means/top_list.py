import scipy.stats as stats
from config import *
from infrastructure.load import *
from infrastructure.save import *
from method.method import *
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn import metrics

num_clusters = 10

method = Method.clustering
clustering_type = ClusteringType.k_means
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
gene_names, gene_vals = get_gene_data(config, gd_type)

clusters_std_avg = []
sil_coeff = []
ch_score = []
for g_id in range(0, len(gene_names)):
    name = gene_names[g_id]
    print(g_id)
    vals = np.array(gene_vals[g_id]).reshape(-1, 1)
    labels = KMeans(n_clusters=num_clusters).fit_predict(vals)
    std_avg = 0.0
    for l_id in range(0, num_clusters):
        cluster_values = []
        for s_id in range(0, len(attributes)):
            if labels[s_id] == l_id:
                cluster_values.append(vals.item(s_id))
        std_avg += np.std(cluster_values)
    std_avg /= float(num_clusters)
    clusters_std_avg.append(std_avg)
    sil_coeff.append(metrics.silhouette_score(vals, labels, metric='euclidean'))
    ch_score.append(metrics.calinski_harabaz_score(vals, labels))

order = np.argsort(clusters_std_avg)
genes_opt = list(np.array(gene_names)[order])
std_avg_opt = list(np.array(clusters_std_avg)[order])
sil_coeff_opt = list(np.array(sil_coeff)[order])
ch_score_opt = list(np.array(ch_score)[order])

fn = 'gene/' + method.value + '/' + clustering_type.value + '_1D_' + gd_type.value + '_' +  geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_features(fn, genes_opt, [std_avg_opt, sil_coeff_opt, ch_score_opt])



