from config.config import *
from infrastructure.load.top import *
import numpy as np

def get_gene_intersection(configs, num_top):
    gene_intersection = []
    if len(configs) > 0:
        gene_intersection = load_top_gene_names(configs[0], num_top)[0:num_top]
        for config in configs[1:]:
            genes = load_top_gene_names(config, num_top)[0:num_top]
            gene_intersection = list(set(gene_intersection).intersection(genes))
    return gene_intersection

def get_gene_intersection_with_article(fn, genes):
    config = Config(read_only=True)
    article_genes = load_top_gene_names_by_article(config, fn)
    intersection_article = list(set(genes).intersection(article_genes))
    return intersection_article