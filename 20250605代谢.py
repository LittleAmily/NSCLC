import os
import gzip
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from seaborn import heatmap
import squidpy as sq
import matplotlib.pyplot as plt
from pathlib import Path
ad = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/stromal.h5ad')
ad  = ad.copy()
ad.obs['stromal_celltype'].cat.categories
cell = ['ADH1B+ CAF', 'Artery', 'BCHE+ SMC', 'COL11A1+ CAF', 'CPE+ Venule',
       'Capillary', 'Lymphatic EC', 'MYH11+ Pericyte', 'Pericyte',
       'SELE+ Venule', 'SMC', 'Tip', 'Venule']
map = {
    'ADH1B+ CAF':'Fibro',
    'Artery':'Artery',
    'BCHE+ SMC':'Fibro',
    'COL11A1+ CAF':'Fibro',
    'CPE+ Venule':'Venule',
    'Capillary':'Capillary',
    'Lymphatic EC':'Lymphatic EC',
    'MYH11+ Pericyte':'Pericyte',
    'SELE+ Venule':'Venule',
    'Pericyte':'Pericyte',
    'SMC':'SMC',
    'Tip':'Tip',
    'Venule':'Venule'
}
ad.obs['stromal_celltype_'] = ad.obs['stromal_celltype'].map(map)
ad.obs['stromal_celltype_'] = ad.obs['stromal_celltype_'].astype('category')
sc.pl.umap(ad, color='stromal_celltype')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype.png', dpi=300, bbox_inches='tight')
sc.tl.leiden(ad,resolution=0.02)
sc.pl.umap(ad, color='leiden')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_leiden.png', dpi=300, bbox_inches='tight')
sc.tl.rank_genes_groups(ad, 'leiden', n_genes=10)
sc.pl.rank_genes_groups(ad, n_genes=10, n_pcs=10, sharey=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/rank_genes_groups.png', dpi=300, bbox_inches='tight')
map = {
    '0':'Artery',
    '1':'Fibro',
    '2':'Pericyte',
    '3':'Venule',
    '4':'Venule',
    '5':'SMC',
    '6':'Lymphatic EC',
    '7':'Capillary',
    '8':'Venule',
    '9':'Fibro',
    '10':'SMC',
    '11':'Fibro'
    
}
ad.obs['stromal_celltype_'] = ad.obs['leiden'].map(map)
ad.obs['stromal_celltype_'] = ad.obs['stromal_celltype_'].astype('category')
sc.pl.umap(ad, color='stromal_celltype_')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype_.png', dpi=1200, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype_.pdf', dpi=300, bbox_inches='tight')
sc.tl.rank_genes_groups(ad, 'stromal_celltype_', n_genes=10)
sc.pl.rank_genes_groups(ad, n_genes=10, n_pcs=10, sharey=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/rank_genes_groups.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/rank_genes_groups.pdf', dpi=300, bbox_inches='tight')

marker_genes = [
                'CD93',     #Artery·
                'PTPRB',     #Artery·
                'VWF',     #Artery·
                'HLA-B',     #B
                'HLA-C',    #CD4
                'ITM2B',    #CD4
                'COL1A2',   #CD8
                'COL6A3',   #CD8
                'COL3A1',    #Treg
                'CCL21',    #Treg
                'TFF3',   #NK
                'TFPI',   #NK
                'CALD1',    #Neutrophil
                'ACTA2',    #Neutrophil
                'TAGLN',   #AM
                'ANK3',   #AM
                'NEDD4L',   #Macro
                'PDE4D',   #Macro
                'SYNE2',    #Mast
                'MUC16',    #Mast
                'ACKR1'
]
# 执行矩阵可视化
sc.pl.matrixplot(
    ad, 
    marker_genes, 
    'stromal_celltype_', 
    standard_scale='var',
    colorbar_title='column scaled\nexpression'
)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/stromal_celltype_sheatmap-.pdf',dpi = 300, bbox_inches='tight')
