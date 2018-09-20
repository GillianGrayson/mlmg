

def clustering_order(clusters):
    curr_cluster = clusters[0]
    mod_cluster = 0
    clusters_sorted = []
    for cl_id in range(0, len(clusters)):
        if clusters[cl_id] != curr_cluster:
            curr_cluster = clusters[cl_id]
            mod_cluster += 1
        clusters_sorted.append(mod_cluster)
    return  clusters_sorted