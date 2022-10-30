# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 20:44:46 2019

@author: Trung Duc Nguyen
"""


#%% (DONE) IMPORT LIBRARY

import pandas as pd
import numpy as np
import using_now as us


#%% (DONE) CALCULATE GCE (GLOBAL CLEAVAGE EFFICIENCY)

def calculate_gce(path):
    df_ctrl = pd.read_csv(f'raw_count/{path}/rawcount.txt',sep ='\t', index_col = 0)
    df_ctrl = df_ctrl[df_ctrl['ctrl'] >= 20]

    df_ctrl_norm = df_ctrl*1000000/df_ctrl.sum()
    df_gce = pd.DataFrame()
    df_gce['d3g1lg4'] = np.log2(df_ctrl_norm['d3g1lg4'] + 0.1) - np.log2(df_ctrl_norm['ctrl'] + 0.1)
    df_gce['d3g2l'] = np.log2(df_ctrl_norm['d3g2l'] + 0.1) - np.log2(df_ctrl_norm['ctrl'] + 0.1)
    
    # normalize vs wt
    df_gce['d3g1lg4'] = df_gce['d3g1lg4'] - df_gce['d3g1lg4']['Variant_006995']
    df_gce['d3g2l'] = df_gce['d3g2l'] - df_gce['d3g2l']['Variant_006995']

    df_gce.columns = ['dro', 'mp']
    return(df_gce)

def merge_gce_2sets():
    
    df_gce1 = calculate_gce('set2')
    df_gce3 = calculate_gce('set1')
    df_gce = pd.concat([df_gce1, df_gce3]).groupby(level=0).mean()
    df_gce.to_csv('./processed_table/gce.txt', sep ='\t',index = True, header = True)

#%% (DONE) CALCULATE LCE RC (LOCAL CLEAVAGE EFFIENCY, CLEAVAGE ACCURACY)

def process_prdt(path, sample):
    # sample = 'd3g1lg4'
    df_ctrl = pd.read_csv(f'raw_count/{path}/rawcount.txt',sep ='\t', index_col = 0)
    df_ctrl = df_ctrl[df_ctrl['ctrl'] >= 20]
    
    df = pd.read_csv(f'raw_count/{path}/{sample}_rawcount.txt', sep ='\t', index_col = 0)
    df = df[df.index.isin(df_ctrl.index)].fillna(0)
    return(df)

def calculate_rc_lce(path):
    # path = 'set2'
    df_ctrl = pd.read_csv(f'raw_count/{path}/rawcount.txt',sep ='\t', index_col = 0)
    df_ctrl = df_ctrl[df_ctrl['ctrl'] >= 20]
    
    df1 = process_prdt(path, 'd3g1lg4')
    df2 = process_prdt(path, 'd3g2l')

    df1_norm = df1 *1000000 / df1.sum().sum()
    df2_norm = df2 *1000000 / df2.sum().sum()
    df_ctrl_norm = df_ctrl*1000000/df_ctrl.sum()
        
    df1_rc =  df1_norm.div(df1_norm.sum(axis = 1), axis =0)
    df2_rc =  df2_norm.div(df2_norm.sum(axis = 1), axis =0)

    df1_norm = df1_norm[df1_norm.index.isin( df1[df1[[str(i) for i in range(-5,6)]].sum(axis =1) >= 15].index)][[str(i) for i in range(-5,6)]]
    df2_norm = df2_norm[df2_norm.index.isin( df2[df2[[str(i) for i in range(-5,6)]].sum(axis =1) >= 15].index)][[str(i) for i in range(-5,6)]]
    
    df1_rc = df1_rc[df1_rc.index.isin( df1[df1[[str(i) for i in range(-5,6)]].sum(axis =1) >= 15].index)][[str(i) for i in range(-5,6)]]
    df2_rc = df2_rc[df2_rc.index.isin( df2[df2[[str(i) for i in range(-5,6)]].sum(axis =1) >= 15].index)][[str(i) for i in range(-5,6)]]

    df1_lce = np.log2((df1_norm + 0.1).divide(df_ctrl_norm['ctrl'] + 0.1, axis='index')).dropna()
    df2_lce = np.log2((df2_norm + 0.1).divide(df_ctrl_norm['ctrl'] + 0.1, axis='index')).dropna()
    
    df1_lce = df1_lce - df1_lce['0']['Variant_006995']
    df2_lce = df2_lce - df2_lce['0']['Variant_006995']
    
    return(df1_rc, df2_rc, df1_lce, df2_lce)    
    # return(df1_rc, df2_rc)    


def merge_lce_rc():
    # df1: d3g1lg4 : mp
    # df2: d3g2l: dro
    
    df1_rc1, df2_rc1, df1_lce1, df2_lce1 = calculate_rc_lce('set2')
    df1_rc3, df2_rc3, df1_lce3, df2_lce3 = calculate_rc_lce('set1')
    
    df1_rc = pd.concat([df1_rc1, df1_rc3]).groupby(level=0).mean()

    df2_rc = pd.concat([df2_rc1, df2_rc3]).groupby(level=0).mean()

    df1_rc.to_csv('./processed_table/mp_rc.txt', sep ='\t',index = True, header = True)
    df2_rc.to_csv('./processed_table/dro_rc.txt', sep ='\t',index = True, header = True)
    
    df1_lce = pd.concat([df1_lce1, df1_lce3]).groupby(level=0).mean()
    df2_lce = pd.concat([df2_lce1, df2_lce3]).groupby(level=0).mean()
    
    df1_lce.to_csv('./processed_table/mp_lce.txt', sep ='\t',index = True, header = True)
    df2_lce.to_csv('./processed_table/dro_lce.txt', sep ='\t',index = True, header = True)

#%% (DONE) CALCULATE HOMO (HOMOGENEITY)

def make_major_cleavage_site():
    
    df1_rc = pd.read_csv('./processed_table/mp_rc.txt', sep ='\t', index_col = 0)
    df2_rc = pd.read_csv('./processed_table/dro_rc.txt', sep ='\t', index_col = 0)
    
    df_major_cleavage = pd.DataFrame(index = list(set(df1_rc.index.tolist() + df2_rc.index.tolist())))
    
    df_major_cleavage['mcl_mp'] = df1_rc.idxmax(axis = 1)
    df_major_cleavage['homo_mp'] = df1_rc.max(axis = 1)
    
    df_major_cleavage['mcl_dro'] = df2_rc.idxmax(axis = 1)
    df_major_cleavage['homo_dro'] = df2_rc.max(axis = 1)
    df_major_cleavage.to_csv('./processed_table/homo.txt', sep ='\t',index = True, header = True)
    
    # df_major_cleavage['second largest homogeneity']  = np.sort(df_rc.values)[:,-2]
    df2_rc.loc['Variant_004622']
    df_major_cleavage.loc['Variant_004622']

