import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import anndata as ad

# 读取数据
data_files = {
    'PA001': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957590_NSCLC_PA001_sn_raw_feature_bc_matrix.h5',
    'PA004': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957591_NSCLC_PA004_sn_raw_feature_bc_matrix.h5',
    'PA005': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957592_NSCLC_PA005_sn_raw_feature_bc_matrix.h5',
    'PA019': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957593_NSCLC_PA019_sn_raw_feature_bc_matrix.h5',
    'PA025': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957594_NSCLC_PA025_sn_raw_feature_bc_matrix.h5',
    'PA034': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957595_NSCLC_PA034_sn_raw_feature_bc_matrix.h5',
    'PA042': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957596_NSCLC_PA042_sn_raw_feature_bc_matrix.h5',
    'PA043': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957597_NSCLC_PA043_sn_raw_feature_bc_matrix.h5',
    'PA048': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957598_NSCLC_PA048_sn_raw_feature_bc_matrix.h5',
    'PA054': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957599_NSCLC_PA054_sn_raw_feature_bc_matrix.h5',
    'PA056': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957600_NSCLC_PA056_sn_raw_feature_bc_matrix.h5',
    'PA060': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957601_NSCLC_PA060_sn_raw_feature_bc_matrix.h5',
    'PA067': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957602_NSCLC_PA067_sn_raw_feature_bc_matrix.h5',
    'PA070': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957604_NSCLC_PA070_sn_raw_feature_bc_matrix.h5',
    'PA072': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957605_NSCLC_PA072_sn_raw_feature_bc_matrix.h5',
    'PA076': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957606_NSCLC_PA076_sn_raw_feature_bc_matrix.h5',
    'PA080': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957607_NSCLC_PA080_sn_raw_feature_bc_matrix.h5',
    'PA104': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957608_NSCLC_PA104_sn_raw_feature_bc_matrix.h5',
    'PA125': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957609_NSCLC_PA125_sn_raw_feature_bc_matrix.h5',
    'PA141': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957610_NSCLC_PA141_sn_raw_feature_bc_matrix.h5',
    'N524': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957611_NSCLC_N254_sn_raw_feature_bc_matrix.h5',
    'N586': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957612_NSCLC_N586_sn_raw_feature_bc_matrix.h5',
    'STK_1': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957613_NSCLC_STK_1_sn_raw_feature_bc_matrix.h5',
    'STK_3': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957614_NSCLC_STK_3_sn_raw_feature_bc_matrix.h5',
    'STK_5dot1': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957615_NSCLC_STK_5dot1_sn_raw_feature_bc_matrix.h5',
    'STK_5dot2': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957616_NSCLC_STK_5dot2_sn_raw_feature_bc_matrix.h5',
    'STK_14': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957617_NSCLC_STK_14_sn_raw_feature_bc_matrix.h5',
    'STK_18': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957618_NSCLC_STK_18_sn_raw_feature_bc_matrix.h5',
    'STK_2': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957619_NSCLC_STK_2_sn_raw_feature_bc_matrix.h5',
    'STK_15': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957620_NSCLC_STK_15_sn_raw_feature_bc_matrix.h5',
    'STK_20': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957621_NSCLC_STK_20_sn_raw_feature_bc_matrix.h5',
    'STK_21': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957622_NSCLC_STK_21_sn_raw_feature_bc_matrix.h5',
    'STK_22dot2': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957623_NSCLC_STK_22dot2_sn_raw_feature_bc_matrix.h5',
    'KRAS_6': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957625_NSCLC_KRAS_6_sn_raw_feature_bc_matrix.h5',
    'KRAS_7': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957626_NSCLC_KRAS_7_sn_raw_feature_bc_matrix.h5',
    'KRAS_8': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957627_NSCLC_KRAS_8_sn_raw_feature_bc_matrix.h5',
    'KRAS_17': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957628_NSCLC_KRAS_17_sn_raw_feature_bc_matrix.h5',
    'KRAS_10': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957629_NSCLC_KRAS_10_sn_raw_feature_bc_matrix.h5',
    'KRAS_11': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957630_NSCLC_KRAS_11_sn_raw_feature_bc_matrix.h5',
    'KRAS_12': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957631_NSCLC_KRAS_12_sn_raw_feature_bc_matrix.h5',
    'KRAS_13': '/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957632_NSCLC_KRAS_13_sn_raw_feature_bc_matrix.h5'
}


PA001 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957590_NSCLC_PA001_sn_raw_feature_bc_matrix.h5')
PA004 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957591_NSCLC_PA004_sn_raw_feature_bc_matrix.h5')
PA005 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957592_NSCLC_PA005_sn_raw_feature_bc_matrix.h5')  
PA019 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957593_NSCLC_PA019_sn_raw_feature_bc_matrix.h5')  
PA025 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957594_NSCLC_PA025_sn_raw_feature_bc_matrix.h5')
PA034 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957595_NSCLC_PA034_sn_raw_feature_bc_matrix.h5')
PA042 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957596_NSCLC_PA042_sn_raw_feature_bc_matrix.h5')
PA043 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957597_NSCLC_PA043_sn_raw_feature_bc_matrix.h5')
PA048 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957598_NSCLC_PA048_sn_raw_feature_bc_matrix.h5')
PA054 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957599_NSCLC_PA054_sn_raw_feature_bc_matrix.h5')
PA056 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957600_NSCLC_PA056_sn_raw_feature_bc_matrix.h5')
PA060 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957601_NSCLC_PA060_sn_raw_feature_bc_matrix.h5')
PA067 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957602_NSCLC_PA067_sn_raw_feature_bc_matrix.h5')
PA070 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957604_NSCLC_PA070_sn_raw_feature_bc_matrix.h5')
PA072 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957605_NSCLC_PA072_sn_raw_feature_bc_matrix.h5')
PA076 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957606_NSCLC_PA076_sn_raw_feature_bc_matrix.h5')
PA080 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957607_NSCLC_PA080_sn_raw_feature_bc_matrix.h5')
PA104 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957608_NSCLC_PA104_sn_raw_feature_bc_matrix.h5')
PA125 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957609_NSCLC_PA125_sn_raw_feature_bc_matrix.h5')
PA141 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957610_NSCLC_PA141_sn_raw_feature_bc_matrix.h5')
N524 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957611_NSCLC_N254_sn_raw_feature_bc_matrix.h5')
N586 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957612_NSCLC_N586_sn_raw_feature_bc_matrix.h5')
STK_1 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957613_NSCLC_STK_1_sn_raw_feature_bc_matrix.h5')
STK_3 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957614_NSCLC_STK_3_sn_raw_feature_bc_matrix.h5')
STK_5dot1 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957615_NSCLC_STK_5dot1_sn_raw_feature_bc_matrix.h5')
STK_5dot2 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957616_NSCLC_STK_5dot2_sn_raw_feature_bc_matrix.h5')
STK_14 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957617_NSCLC_STK_14_sn_raw_feature_bc_matrix.h5')
STK_18 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957618_NSCLC_STK_18_sn_raw_feature_bc_matrix.h5')
STK_2 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957619_NSCLC_STK_2_sn_raw_feature_bc_matrix.h5')
STK_15 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957620_NSCLC_STK_15_sn_raw_feature_bc_matrix.h5')
STK_20 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957621_NSCLC_STK_20_sn_raw_feature_bc_matrix.h5')
STK_21 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957622_NSCLC_STK_21_sn_raw_feature_bc_matrix.h5')
STK_22dot2 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957623_NSCLC_STK_22dot2_sn_raw_feature_bc_matrix.h5')
KRAS_6 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957625_NSCLC_KRAS_6_sn_raw_feature_bc_matrix.h5')
KRAS_7 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957626_NSCLC_KRAS_7_sn_raw_feature_bc_matrix.h5')
KRAS_8 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957627_NSCLC_KRAS_8_sn_raw_feature_bc_matrix.h5')
KRAS_17 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957628_NSCLC_KRAS_17_sn_raw_feature_bc_matrix.h5')
KRAS_10 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957629_NSCLC_KRAS_10_sn_raw_feature_bc_matrix.h5')
KRAS_11 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957630_NSCLC_KRAS_11_sn_raw_feature_bc_matrix.h5')
KRAS_12 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957631_NSCLC_KRAS_12_sn_raw_feature_bc_matrix.h5')
KRAS_13 = sc.read_10x_h5('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499/GSM6957632_NSCLC_KRAS_13_sn_raw_feature_bc_matrix.h5')

PA001.obs_names = ['PA001_' + x for x in PA001.obs_names]
PA004.obs_names = ['PA004_' + x for x in PA004.obs_names]
PA005.obs_names = ['PA005_' + x for x in PA005.obs_names]
PA019.obs_names = ['PA019_' + x for x in PA019.obs_names]
PA025.obs_names = ['PA025_' + x for x in PA025.obs_names]
PA034.obs_names = ['PA034_' + x for x in PA034.obs_names]
PA042.obs_names = ['PA042_' + x for x in PA042.obs_names]
PA043.obs_names = ['PA043_' + x for x in PA043.obs_names]
PA048.obs_names = ['PA048_' + x for x in PA048.obs_names]
PA054.obs_names = ['PA054_' + x for x in PA054.obs_names]
PA056.obs_names = ['PA056_' + x for x in PA056.obs_names]
PA060.obs_names = ['PA060_' + x for x in PA060.obs_names]
PA067.obs_names = ['PA067_' + x for x in PA067.obs_names]
PA070.obs_names = ['PA070_' + x for x in PA070.obs_names]
PA072.obs_names = ['PA072_' + x for x in PA072.obs_names]
PA076.obs_names = ['PA076_' + x for x in PA076.obs_names]
PA080.obs_names = ['PA080_' + x for x in PA080.obs_names]
PA104.obs_names = ['PA104_' + x for x in PA104.obs_names]
PA125.obs_names = ['PA125_' + x for x in PA125.obs_names]
PA141.obs_names = ['PA141_' + x for x in PA141.obs_names]
N524.obs_names = ['N524_' + x for x in N524.obs_names]
N586.obs_names = ['N586_' + x for x in N586.obs_names]
STK_1.obs_names = ['STK_1_' + x for x in STK_1.obs_names]
STK_3.obs_names = ['STK_3_' + x for x in STK_3.obs_names]
STK_5dot1.obs_names = ['STK_5dot1_' + x for x in STK_5dot1.obs_names]
STK_5dot2.obs_names = ['STK_5dot2_' + x for x in STK_5dot2.obs_names]
STK_14.obs_names = ['STK_14_' + x for x in STK_14.obs_names]
STK_18.obs_names = ['STK_18_' + x for x in STK_18.obs_names]
STK_2.obs_names = ['STK_2_' + x for x in STK_2.obs_names]
STK_15.obs_names = ['STK_15_' + x for x in STK_15.obs_names]
STK_20.obs_names = ['STK_20_' + x for x in STK_20.obs_names]
STK_21.obs_names = ['STK_21_' + x for x in STK_21.obs_names]
STK_22dot2.obs_names = ['STK_22dot2_' + x for x in STK_22dot2.obs_names]
KRAS_6.obs_names = ['KRAS_6_' + x for x in KRAS_6.obs_names]
KRAS_7.obs_names = ['KRAS_7_' + x for x in KRAS_7.obs_names]
KRAS_8.obs_names = ['KRAS_8_' + x for x in KRAS_8.obs_names]
KRAS_17.obs_names = ['KRAS_17_' + x for x in KRAS_17.obs_names]
KRAS_10.obs_names = ['KRAS_10_' + x for x in KRAS_10.obs_names]
KRAS_11.obs_names = ['KRAS_11_' + x for x in KRAS_11.obs_names]
KRAS_12.obs_names = ['KRAS_12_' + x for x in KRAS_12.obs_names]
KRAS_13.obs_names = ['KRAS_13_' + x for x in KRAS_13.obs_names]
sample_sheet = [PA001, PA004, PA005, PA019, PA025, PA034, PA042, PA043, PA048, PA054, PA056, 
                PA060, PA067, PA070,PA072, PA076, PA080, PA104, PA125, PA141, N524, N586, 
                STK_1, STK_3, STK_5dot1, STK_5dot2, STK_14, STK_18, STK_2, STK_15, STK_20, 
                STK_21, STK_22dot2, KRAS_6, KRAS_7, KRAS_8, KRAS_17, KRAS_10, KRAS_11, KRAS_12, KRAS_13]
for adata in sample_sheet:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
Tagore_S_2025 = ad.concat([PA001, PA004, PA005, PA019, PA025, PA034, PA042, PA043, PA048, PA054, 
                           PA056, PA060, PA067, PA070, PA072, PA076, PA080, PA104, PA125, PA141, 
                           N524, N586, STK_1, STK_3, STK_5dot1, STK_5dot2, STK_14, STK_18, STK_2, 
                           STK_15, STK_20, STK_21, STK_22dot2, KRAS_6, KRAS_7, KRAS_8, KRAS_17, 
                           KRAS_10, KRAS_11, KRAS_12, KRAS_13], label="sample")
integrated_data = pd.read_csv('/home/data/sdzl14/NSCLC/Tagore_S_2025/GSE223499_sn_integrated_data.csv')
# 提取 Unnamed: 0 列并重命名为 UMI
integrated_data['UMI'] = integrated_data['Unnamed: 0']
integrated_data.drop(columns=['Unnamed: 0'], inplace=True)
QC_UMI = integrated_data['UMI']
QC_UMI = QC_UMI.tolist()
Tagore_S_2025 = Tagore_S_2025[Tagore_S_2025.obs_names.isin(QC_UMI)]
common_umi = list(set(Tagore_S_2025.obs_names).intersection(set(integrated_data['UMI'])))
integrated_data = integrated_data[integrated_data['UMI'].isin(common_umi)]
integrated_data.set_index('UMI', inplace=True)
Tagore_S_2025 = Tagore_S_2025[Tagore_S_2025.obs_names.isin(common_umi)]

# 将 filtered_integrated_data 中的列添加到 Tagore_S_2025.obs 中
# 将 filtered_integrated_data 中的列添加到 Tagore_S_2025.obs 中
Tagore_S_2025.obs = Tagore_S_2025.obs.join(integrated_data)
Tagore_S_2025.write_h5ad('/home/data/sdzl14/NSCLC/Tagore_S_2025/Tagore_S_2025_integrated_data.h5ad')
Peng_Zhang = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/NSCLC-singlecell-all_new.h5ad')

Peng_Zhang = pd.read_csv('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/counts.csv')
meta = pd.read_csv(
   "/home/data/sdzl14/NSCLC/Peng_Zhang_2024/meta.data.csv", index_col="Sample"
)
ad = sc.read_h5ad('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/10X_Lung_Tumour_Annotated_v2.h5ad')
ad.X.max()
sc.pl.umap(ad, color="Cell types")
plt.savefig('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/fig/Cell types.png', dpi=300, bbox_inches='tight')
sc.pl.umap(ad, color="Cell types v25")
plt.savefig('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/fig/Cell types v25.png', dpi=300, bbox_inches='tight')
sc.pl.umap(ad, color="leiden")
plt.savefig('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/fig/leiden.png', dpi=300, bbox_inches='tight')
sc.pl.umap(ad, color="PHASE")
plt.savefig('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/fig/PHASE.png', dpi=300, bbox_inches='tight')
ad = sc.read_h5ad('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/10X_Lung_Healthy_Background_Annotated_v2.h5ad')
adata = sc.AnnData(Peng_Zhang.T)  # Scanpy期望基因在列上，所以转置矩阵
adata[1,]
meta.drop(columns=['Unnamed: 0'], inplace=True)
adata.write_h5ad('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/scRNA_annotated.h5ad')
print(Peng_Zhang.X[:5, :5]) 

from scipy.sparse import csr_matrix
adata.X = csr_matrix(adata.X)
print(type(Peng_Zhang.X))
print(adata.X[:5, :5])
adata.write_h5ad('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/scRNA_annotated.h5ad')
ad.var = pd.DataFrame(index=ad.var.index)
ad.obsm = {}
ad.write_h5ad('/home/data/sdzl14/NSCLC/EBI_De_Zuani_2024/SC-data/10X_Lung_Healthy_Background_Annotated.h5ad')
sc111 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/1.1.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)
sc121 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/1.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequentreading
)
sc122 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/1.2.2/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for fastersubsequent reading
)
sc131 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/1.3.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfaster subsequentreading
)
sc141 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/1.4.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequentreading
)
sc211 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/2.1.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequentreading
)
sc221 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/2.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfaster subsequentreading
)
sc222 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/2.2.2/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequentreading
)
sc231 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfaster subsequentreading
)
sc241 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfaster subsequentreading
)
sc251 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/5.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc311 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.1.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc321 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc322 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.2.2/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc323 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.2.3/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc331 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.3.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc341 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.4.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc351 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/3.5.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc411 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.1.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc421 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc422 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.2.2/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc431 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.3.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc441 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.4.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc451 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/4.5.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc521 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/5.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc522 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/5.2.2/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc531 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/5.3.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc551 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/5.5.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc621 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/6.2.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc631 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/6.3.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading    
)
sc641 = sc.read_10x_mtx(
    "/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_data_count/6.4.1/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file forfastersubsequentreading
)
sc111.obs_names = ['1.1.1_' + x for x in sc111.obs_names]
sc121.obs_names = ['1.2.1_' + x for x in sc121.obs_names]
sc122.obs_names = ['1.2.2_' + x for x in sc122.obs_names]
sc131.obs_names = ['1.3.1_' + x for x in sc131.obs_names]
sc141.obs_names = ['1.4.1_' + x for x in sc141.obs_names]
sc211.obs_names = ['2.1.1_' + x for x in sc211.obs_names]
sc221.obs_names = ['2.2.1_' + x for x in sc221.obs_names]
sc222.obs_names = ['2.2.2_' + x for x in sc222.obs_names]
sc231.obs_names = ['2.3.1_' + x for x in sc231.obs_names]
sc241.obs_names = ['2.4.1_' + x for x in sc241.obs_names]
sc251.obs_names = ['2.5.1_' + x for x in sc251.obs_names]
sc311.obs_names = ['3.1.1_' + x for x in sc311.obs_names]
sc321.obs_names = ['3.2.1_' + x for x in sc321.obs_names]
sc322.obs_names = ['3.2.2_' + x for x in sc322.obs_names]
sc323.obs_names = ['3.2.3_' + x for x in sc323.obs_names]
sc331.obs_names = ['3.3.1_' + x for x in sc331.obs_names]
sc341.obs_names = ['3.4.1_' + x for x in sc341.obs_names]
sc351.obs_names = ['3.5.1_' + x for x in sc351.obs_names]
sc411.obs_names = ['4.1.1_' + x for x in sc411.obs_names]
sc421.obs_names = ['4.2.1_' + x for x in sc421.obs_names]
sc422.obs_names = ['4.2.2_' + x for x in sc422.obs_names]
sc431.obs_names = ['4.3.1_' + x for x in sc431.obs_names]
sc441.obs_names = ['4.4.1_' + x for x in sc441.obs_names]
sc451.obs_names = ['4.5.1_' + x for x in sc451.obs_names]
sc521.obs_names = ['5.2.1_' + x for x in sc521.obs_names]
sc522.obs_names = ['5.2.2_' + x for x in sc522.obs_names]
sc531.obs_names = ['5.3.1_' + x for x in sc531.obs_names]
sc551.obs_names = ['5.5.1_' + x for x in sc551.obs_names]
sc621.obs_names = ['6.2.1_' + x for x in sc621.obs_names]
sc631.obs_names = ['6.3.1_' + x for x in sc631.obs_names]
sc641.obs_names = ['6.4.1_' + x for x in sc641.obs_names]

sample_sheet = [
    sc111, sc121, sc122, sc131, sc141,
    sc211, sc221, sc222, sc231, sc241, sc251,
    sc311, sc321, sc322, sc323, sc331, sc341, sc351,
    sc411, sc421, sc422, sc431, sc441, sc451,
    sc521, sc522, sc531, sc551,
    sc621, sc631, sc641
]
for adata in sample_sheet:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
Peilin_2025 = ad.concat(sample_sheet, label="sample")
meta = pd.read_csv('~/NSCLC/Peilin_Wang_2025/meta.data.csv')
# 提取 Unnamed: 0 列并重命名为 UMI
meta['UMI'] = meta['Unnamed: 0']
meta.drop(columns=['Unnamed: 0'], inplace=True)
QC_UMI = meta['UMI']
QC_UMI = QC_UMI.tolist()
Peilin_2025 = Peilin_2025[Peilin_2025.obs_names.isin(QC_UMI)]
common_umi = list(set(Peilin_2025.obs_names).intersection(set(meta['UMI'])))
meta = meta[meta['UMI'].isin(common_umi)]
meta.set_index('UMI', inplace=True)
Peilin_2025 = Peilin_2025[Peilin_2025.obs_names.isin(common_umi)]

# 将 filtered_integrated_data 中的列添加到 Tagore_S_2025.obs 中
# 将 filtered_integrated_data 中的列添加到 Tagore_S_2025.obs 中
Peilin_2025.obs = Peilin_2025.obs.join(meta)
Peilin_2025.write_h5ad('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_annotated.h5ad')

MSK1263_A1 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_A1_RNA')
MSK1263_A3 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_A3_RNA')
MSK1263_A4 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_A4_RNA')  
MSK1263_LN = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_LN_RNA')
MSK1263_Normal = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_Normal_RNA')
MSK1263_R1 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R1_RNA')
MSK1263_R2 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R2_RNA')
MSK1263_R3 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R3_RNA')
MSK1263_R4 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R4_RNA')
MSK1263_R5 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R5_RNA')
MSK1263_R6 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R6_RNA')
MSK1263_R7 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R7_RNA')
MSK1263_R8 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1263_R8_RNA')
MSK1302_LN = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_LN_RNA')
MSK1302_Normal = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_Normal_RNA')
MSK1302_R1 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R1_RNA')
MSK1302_R2 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R2_RNA')
MSK1302_R3 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R3_RNA')
MSK1302_R4 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R4_RNA')
MSK1302_R5 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R5_RNA')
MSK1302_R6 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R6_RNA')
MSK1302_R7 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R7_RNA')
MSK1302_R8 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1302_R8_RNA')
MSK1344_Normal = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_Normal_RNA')
MSK1344_R1 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R1_RNA')
MSK1344_R2 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R2_RNA')
MSK1344_R3 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R3_RNA')
MSK1344_R4 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R4_RNA')
MSK1344_R5 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R5_RNA')
MSK1344_R6 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R6_RNA')
MSK1344_R7 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R7_RNA')
MSK1344_R8 = sc.read_10x_mtx('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/10_matrix/MSK1344_R8_RNA')
sample_sheet = [
    MSK1263_A1, MSK1263_A3, MSK1263_A4, MSK1263_LN, MSK1263_Normal,
    MSK1263_R1, MSK1263_R2, MSK1263_R3, MSK1263_R4, MSK1263_R5,
    MSK1263_R6, MSK1263_R7, MSK1263_R8,
    MSK1302_LN, MSK1302_Normal, MSK1302_R1, MSK1302_R2, MSK1302_R3,
    MSK1302_R4, MSK1302_R5, MSK1302_R6, MSK1302_R7, MSK1302_R8,
    MSK1344_Normal, MSK1344_R1, MSK1344_R2, MSK1344_R3, MSK1344_R4,
    MSK1344_R5, MSK1344_R6, MSK1344_R7, MSK1344_R8
]
for adata in sample_sheet:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
# 定义样本名称列表（需与 sample_sheet 中对象顺序一致）
sample_names = [
    'MSK1263_A1', 'MSK1263_A3', 'MSK1263_A4', 'MSK1263_LN', 'MSK1263_Normal',
    'MSK1263_R1', 'MSK1263_R2', 'MSK1263_R3', 'MSK1263_R4', 'MSK1263_R5',
    'MSK1263_R6', 'MSK1263_R7', 'MSK1263_R8',
    'MSK1302_LN', 'MSK1302_Normal', 'MSK1302_R1', 'MSK1302_R2', 'MSK1302_R3',
    'MSK1302_R4', 'MSK1302_R5', 'MSK1302_R6', 'MSK1302_R7', 'MSK1302_R8',
    'MSK1344_Normal', 'MSK1344_R1', 'MSK1344_R2', 'MSK1344_R3', 'MSK1344_R4',
    'MSK1344_R5', 'MSK1344_R6', 'MSK1344_R7', 'MSK1344_R8'
]

# 执行批量处理
for name, adata in zip(sample_names, sample_sheet):
    adata.obs_names = [f"{name}_{x}" for x in adata.obs_names]
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
import scrublet as scr
for adata in sample_sheet:
    
    ad = adata.copy()
    scrub = scr.Scrublet(ad.X, expected_doublet_rate=ad.n_obs*8*1e-6)
    out = scrub.scrub_doublets(verbose=False, n_prin_comps=20)
    scrdf = pd.DataFrame({'scr_score': out[0], 'scr_pred': out[1]}, index=ad.obs.index)
    ad.obs = pd.concat([ad.obs, scrdf], axis=1)
    adata = ad.copy()
import anndata as ad    
Ansuman_Satpathy_2023 = ad.concat(sample_sheet, label='sample')
Ansuman_Satpathy_2023.obs_names



Ansuman_Satpathy_2023.write_h5ad('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023.h5ad')
Ansuman_Satpathy_2023 = sc.read('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023.h5ad')
# 识别线粒体、核糖体和红血球基因
# mitochondrial genes
Ansuman_Satpathy_2023.var['mt'] = Ansuman_Satpathy_2023.var_names.str.startswith('MT-')
# ribosomal genes
Ansuman_Satpathy_2023.var['ribo'] = Ansuman_Satpathy_2023.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
Ansuman_Satpathy_2023.var['hb'] = Ansuman_Satpathy_2023.var_names.str.contains(('HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ'))

# 计算质控指标
sc.pp.calculate_qc_metrics(Ansuman_Satpathy_2023, qc_vars=['mt', 'ribo', 'hb'],
                           percent_top=None, log1p=False, inplace=True)
Ansuman_Satpathy_2023.write_h5ad('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023_raw.h5ad', compression='gzip')

# 创建数字到名称的映射字典（0-31对应32个样本）
sample_mapping = {str(i): name for i, name in enumerate(sample_names)}  # 关键修改：键转为字符串

# 执行映射转换
Ansuman_Satpathy_2023.obs['sample'] = Ansuman_Satpathy_2023.obs['sample'].map(sample_mapping)







### 查看过滤前的质控指标
sc.pl.violin(Ansuman_Satpathy_2023, ['n_genes_by_counts', 'total_counts'], jitter=0.4,
             groupby = 'sample', rotation= 45)
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/QC-1.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/QC-1.png',dpi = 300, bbox_inches='tight')
sc.pl.violin(Ansuman_Satpathy_2023, ['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'], jitter=0.4,
             groupby = 'sample', rotation= 45)
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/QC-2.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/QC-2.png',dpi = 300, bbox_inches='tight')
sc.pp.filter_cells(Ansuman_Satpathy_2023, min_genes=450)
sc.pp.filter_genes(Ansuman_Satpathy_2023, min_cells=1)
#Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['n_genes_by_counts'] <= 4000, :]
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['total_counts'] <= 15000, :]
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['pct_counts_mt'] <= 15, :]
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['pct_counts_hb'] <= 1, :]
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['sample'] != 'MSK1344_R1']
Ansuman_Satpathy_2023
sc.pl.violin(Ansuman_Satpathy_2023, ['n_genes_by_counts', 'total_counts'], jitter=0.4,
             groupby = 'sample', rotation= 45)
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/after_QC-1.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/after_QC-1.png',dpi = 300, bbox_inches='tight')

sc.pl.violin(Ansuman_Satpathy_2023, ['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'], jitter=0.4,
             groupby = 'sample', rotation= 45)
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/after_QC_p2.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/QC/after_QC_p2.png',dpi = 300, bbox_inches='tight')


# 5. 保存过滤后的数据
Ansuman_Satpathy_2023.write_h5ad(
    '/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023_filtered.h5ad'
)
########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 1. 数据整合与预处理（使用Scanpy实现类似Seurat流程）
# 排除TCR基因（假设TCR相关基因为TRAC, TRBC, TRDC等）
Ansuman_Satpathy_2023 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023_filtered.h5ad')
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023.copy()
tcr_genes = ['TRAC', 'TRBC', 'TRDC', 'TRGV', 'TRDV', 'TRDV']
Ansuman_Satpathy_2023 = Ansuman_Satpathy_2023[:, ~Ansuman_Satpathy_2023.var_names.isin(tcr_genes)]
# 校正测序深度
sc.pp.normalize_total(Ansuman_Satpathy_2023)
sc.pp.log1p(Ansuman_Satpathy_2023)
sc.pp.highly_variable_genes(Ansuman_Satpathy_2023, n_top_genes=2000, flavor='seurat')

import scanorama
# 按照分组拆分单细胞数据为list
batch_cats = Ansuman_Satpathy_2023.obs['sample'].cat.categories.tolist()
Ansuman_Satpathy_2023_list = [Ansuman_Satpathy_2023[Ansuman_Satpathy_2023.obs['sample'] == b].copy() for b in batch_cats]
# 运行Scanorama对list对象进行整合分析
scanorama.integrate_scanpy(Ansuman_Satpathy_2023_list)

Ansuman_Satpathy_2023.obsm["X_scanorama"] = np.zeros((Ansuman_Satpathy_2023.shape[0], Ansuman_Satpathy_2023_list[0].obsm["X_scanorama"].shape[1]))
for i, b in enumerate(batch_cats):
    Ansuman_Satpathy_2023.obsm["X_scanorama"][Ansuman_Satpathy_2023.obs['sample'] == b] = Ansuman_Satpathy_2023_list[i].obsm["X_scanorama"]
Ansuman_Satpathy_2023



# 2. 维度约简与聚类
sc.pp.neighbors(Ansuman_Satpathy_2023, use_rep="X_scanorama")
sc.tl.umap(Ansuman_Satpathy_2023)

sc.pl.embedding(Ansuman_Satpathy_2023,basis="umap", color=['sample'], 
                frameon=False, title='Scanorama',
                legend_loc='right margin', legend_fontsize=12,
                legend_fontweight='normal', legend_fontoutline=6,
                palette=None, cmap=None
)
results_file = '/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Fig/Integrate/'
sc.plt.savefig(results_file + 'umap_scanorama.png', dpi=300, bbox_inches='tight')