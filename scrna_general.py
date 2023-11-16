import os
import re
import scanpy as sc 
import pandas as pd
sc.logging.print_header()
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=81, facecolor='white')

# optional: get prefix of every sample and load to annadata,attention the file name difference between cellranger version


directory = "/data/dk/papillary_thyroid_carcinoma/raw/"

files = os.listdir(directory)
prefixes = []
if any('genes' in item for item in files):
    for file in files:
        match = re.match(r"(.*)genes\.tsv",file)
        if match:
            prefix = match.group(1)
            print(prefix)
            prefixes.append(prefix)
else:
    for file in files:
        match = re.match(r"(.*)features\.tsv",file)
        if match:
            prefix = match.group(1)
            print(prefix)
            prefixes.append(prefix)
            new_filename = file.replace('features', 'genes')
            print(new_filename)
            os.rename(os.path.join(directory, file), os.path.join(directory, new_filename))


adatas = [sc.read_10x_mtx(directory, cache=True, var_names='gene_symbols', prefix=prefix) for prefix in prefixes]
adata = adatas[0].concatenate(adatas[1:])

#preprocess
adata.var_names_make_unique()#deduplicated gene sysmbols removal
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata #freezes anndata state
adata = adata[:, adata.var.highly_variable] 
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

#neighborring and cluster
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=19)
sc.tl.umap(adata)
sc.tl.leiden(adata,resolution=0.1)
sc.pl.umap(adata, color='leiden')

#find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
markers = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
index = markers.values.flatten().tolist()
adata.raw[:,index]