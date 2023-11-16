import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=81,facecolor='white')
adata_sham1 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319987_sham1_')
adata_sham2 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319988_sham2_')
adata_sham3 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319989_sham3_')
adata_MCAO1 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319990_MCAO1_')
adata_MCAO2 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319991_MCAO2_')
adata_MCAO3 = sc.read_10x_mtx('/data/dk/ischematic_stroke/',
                        var_names='gene_symbols',
                        cache=True,
                        prefix='GSM5319992_MCAO3_')


adata_sham1.obs['label'] = 'sham1'
adata_sham2.obs['label'] = 'sham2'
adata_sham3.obs['label'] = 'sham3'
adata_MCAO1.obs['label'] = 'MCAO1'
adata_MCAO2.obs['label'] = 'MCAO2'
adata_MCAO3.obs['label'] = 'MCAO3'

adata_sham1.obs['group'] = 'sham'
adata_sham2.obs['group'] = 'sham'
adata_sham3.obs['group'] = 'sham'
adata_MCAO1.obs['group'] = 'MCAO'
adata_MCAO2.obs['group'] = 'MCAO'
adata_MCAO3.obs['group'] = 'MCAO'
adata = adata_sham1.concatenate(adata_sham1,adata_sham2,adata_sham3,adata_MCAO1,adata_MCAO2,adata_MCAO3)
adata.var_names_make_unique() 
adata
#preprocess data,filter contamination data
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  #
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
#normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#focus on high variable genes 
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
#raw data
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
#pca 
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

#computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=19)
#embedding the neighborhood graph
sc.tl.paga(adata,groups='label')
sc.pl.paga(adata)
sc.tl.umap(adata,init_pos='paga')
sc.pl.umap(adata)
sc.tl.leiden(adata,resolution=0.5)
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='group')
#find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
marker = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
marker.to_csv('/data/dk/ischematic_stroke/result/marker.csv',index=False)
marker.iloc[:,1]
#','.join(marker.iloc[:,24].astype(str))
new_cluster_names = ['Foam cell','Pyramidal neuron','Microglia cell','CD14+ monocyte','Vascular cell','Choroid cell','Quiescent neural stem cell','Mural cell','Quiescent neural stem cell2','Antigen-presenting cell','Mature macrophage','Schwann cell','erythrocyte','Pericyte-like cell','Quiescent neural stem cell3','Monocyte','Microglia cell2','Lymphocyte','Leptomeningeal cell','Hippocampal pyramidal precursor cell','Ependymal cell','Activated neural stem cell','Neuroblast','HSPC']
adata.rename_categories('leiden', new_cluster_names)
#plot umap of annotated cell cluster
sc.pl.umap(adata, color='leiden', legend_loc='right margin',legend_fontsize=6, title='', frameon=False)

#get cell with the largest numbers of differential genes in mcao group

#get diff genes of mcao vs sham
sc.tl.rank_genes_groups(adata,groupby='group',groups=['MCAO','sham'], 
                        reference='sham',method='wilcoxon', key_added='condition')
diff_genes = sc.get.rank_genes_groups_df(adata,group='MCAO',pval_cutoff=0.01,key='condition',log2fc_min=0)['names'].values
#get diff gene of every cell line, and get the count numbers of 
#intersect with different genes between maco and sham.
len(diff_genes)
df = pd.DataFrame()
for fields in adata.uns['rank_genes_groups']['names'].dtype.fields.keys():
    print(fields)
    col = sc.get.rank_genes_groups_df(adata,group=fields,pval_cutoff=0.01,log2fc_min=1)['names'].values
    print(len(np.intersect1d(diff_genes,col)))
    if fields == 'Ependymal cell':
        genelist = np.intersect1d(diff_genes,col)
len(genelist)
# enrcihment analysis 
import gseapy as gp
glist = genelist.tolist()
print(glist)
enr = gp.enrichr(gene_list=glist,
                 gene_sets=['KEGG_2019_Mouse'],
                 organism='mouse',)


enr.results.head(5)
from gseapy import barplot, dotplot
ax = dotplot(enr.results,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=10,
              top_term=5,
              figsize=(3,5),
              title = "KEGG",
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
             )




import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata,color='leiden',legend_loc='on data')
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata,color='leiden',legend_loc='right margin')






