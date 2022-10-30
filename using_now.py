# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:58:24 2021

@author: Trung Duc Nguyen
"""

#%% IMPORT LIBRARY

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools
import common_tools as ct
import re
import statsmodels.api as sm

#%% (DONE) TOOLs

def make_motif(pos1, pos2, pos3, pos4, pos5, pos6):
    motif = ''
    for pos in [pos1, pos2, pos3, pos4, pos5, pos6]:
        if pos == 'G-T':
            motif = motif +'W'
        elif pos == 'T-G':
            motif = motif +'w'
        elif pos in ['A-T', 'T-A', 'G-C', 'C-G']:
            motif = motif +pos[0]
        elif pos in ['A-C', 'C-A']:
            motif = motif +pos[0].lower()
        else:
            motif = motif + 'N'
    return(motif)


def call_mother_of_motif(motif, number_nt):
    # motif = 'pos1 C-G; pos2 A-T; pos3 T-A; pos4 G-C; pos5 T-A; pos6 A-T'
    # number_nt  = 4
    
    list_motif = [i for i in motif[:-2].split('; ')]
    
    list_mother = []
        
    for mother in itertools.combinations(list_motif, number_nt):
        # print('; '.join(list(mother)))
        list_mother.append('; '.join(list(mother)) + '; ')
    return(list_mother)

def extract_variant_given_sg (df_sg, list_sg):
    list_variant = []
    for sg in list_sg:
        list_variant = list_variant + df_sg[df_sg[sg] != 'no'].index.tolist()
    return(list(set(list_variant)))

def extract_variant_given_motif(df_input, motif):
    df_out = df_input.copy() 
    for pos_bp in motif.split('; '):
        if pos_bp != '':
            pos = pos_bp.split(' ')[0]
            bp = pos_bp.split(' ')[1]
            df_out = df_out[df_out[pos] == bp]
    return(df_out)


def extract_list_symetric_struct(df_all):
    # mismatch included
    
    df_sub_struct = df_all.copy()
    
    df_sub_struct = df_sub_struct[df_sub_struct['mian'].str.contains('1-FFFFFFFFFFFFFFFFFF-18')]
    df_sub_struct = df_sub_struct[df_sub_struct['struct'].str.contains('91-UUUUUUUU-98 99-MMLLLLMM-106 107-TT-108')]
    df_sub_struct = df_sub_struct[~df_sub_struct['struct'].str.contains('A')]
    df_sub_struct = df_sub_struct[~df_sub_struct['struct'].str.contains('B')]
    list_symetric_struct = df_sub_struct['name_struct'].unique().tolist()

    return(list_symetric_struct)

def draw_heatmap_cleavage_accuracy(df_draw, save_fig):
    # df_draw = df_rc_mp.copy()
    # df_draw = df_rc_dro.copy()
    
    plt.figure(figsize = (4,1))
    ax = sns.heatmap(df_draw.sort_values(by= ['0', '-1', '1'], ascending=True), 
                 yticklabels=False, xticklabels=True, #columns, 
                 cmap = sns.color_palette('GnBu', 200), vmax = 0.5, vmin = 0, cbar=False,)
    
    ct.decorate_boxplot(ax)

    for i in range(1, 11):
        plt.axvline(i, c="white", linewidth=2)
    
    plt.ylabel('', size = 0)
    plt.xlabel('', size = 0)
    
    plt.xticks( color ='white', size =0)
    plt.xticks([])
    plt.yticks([])

    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
        spine.set_color('black')

    plt.savefig(f'{save_fig} {len(df_draw)} variants' , bbox_inches="tight", dpi =300)


def draw_heatmap_rc_given_variant(list_variant, df_rc, save_fig):
    draw_heatmap_cleavage_accuracy(df_rc[df_rc.index.isin(list_variant)], save_fig)



def draw_heatmap_rc_pattern_given_pos_bp(df_using, df_rc, list_pos, list_bp, save_fig, draw = False):
    df = df_using.copy()
    for pos, bp in zip(list_pos, list_bp):
        pos  = 'pos' + str(pos)
        df = df[df[pos].isin(bp)]
        save_fig = save_fig + '{} {}; '.format(pos, ' '.join(bp))
        
        
    df_draw = df_rc[df_rc.index.isin(df.index)]
    # df_draw = df_rc14[df_rc14.index.isin(df_using.dropna().index)]


    if draw == True:
        draw_heatmap_cleavage_accuracy(df_draw, save_fig )
    # draw_heatmap_cleavage_accuracy(df_draw, save_fig + str(len(df_draw)) + ' d3g1lg4')
    return(df_draw)

def make_new_sequence_code_based_on_struct(seq, struct):
    # seq = 'CTGCCATTTTACAATCCAACAAAAAATCTAATTTCTCCACGTCTTTGGTAATAAGGTTTGGCAAAGATGTGGAAAAATTGGACTTTTCGTCTTTCAACCCAACCGGAA'
    # struct = 'FFFFFFFFFFFFFFFFFFMMSMMMMSMMwMMMMMSMMMMMWMMMMMMSwMMLLLLLLMMWSMMMMMMwMMMMMSMMMMMWMMSMMMMSMMUUUUUUUUMMLLLLMMTT'
    
    new_seq = ''
    for i in range(len(seq)):
        if struct[i] in ['F', 'S', 'A', 'U', 'T', 'B']:
            new_seq = new_seq + '.'
        elif struct[i] in ['L']:
            new_seq = new_seq + '-'
        elif struct[i] in ['w', 'W']:
            new_seq = new_seq + struct[i]
        else:
            new_seq = new_seq + seq[i]
    return(new_seq)

#%% (DONE) IMPORT STRUCTURE,SG, SEQUENCE

def import_refseq():
    df_refseq = pd.read_csv('./processed_table/576_ref.txt', sep ='\t', names = ['var', 'ref_seq']).set_index('var')
    return(df_refseq)

def import_sg():
    df_sg = pd.read_csv('processed_table/576_ref_allSG_2nd.txt', sep ='\t').set_index('var')
    return(df_sg)

def import_struct():
    df_struct = pd.read_csv('processed_table/576_structure_concrete.txt', sep ='\t').set_index('Variant')
    
    df_name = pd.DataFrame()
    df_name['struct'] = df_struct['New_define_structure'].value_counts()
    df_name['Name_struct'] = ['struct_' + str(i).zfill(4) for i in range(1, 3012)]
    df_struct = pd.merge(df_struct, df_name[['Name_struct']], left_on ='New_define_structure', right_index = True)
    return(df_struct)

def import_struct_with_classification():
    df_struct = pd.read_csv('processed_table/576_struct_with_classification.txt', sep ='\t', index_col = 0).fillna('')
    return(df_struct)


#%% (DONE) COMBINE COUNT, SCORE INTO DF_ALL

def import_ctrl_count(path):
    df_ctrl = pd.read_csv(f'./raw_count/{path}/ctrl_sumcount.txt',
                          sep ='\t', names = ['var', 'barcode_number', 'count'] ).set_index('var')
    return(df_ctrl)

def import_raw_count(path):
    df_ctrl = pd.read_csv(f'./raw_count/{path}/rawcount.txt',
                          sep ='\t').set_index('var')
    return(df_ctrl)

def import_gce():
    df_gce = pd.read_csv('processed_table/gce.txt', sep ='\t').set_index('var')
    return(df_gce)

def import_rc(sample):
    df_rc = pd.read_csv('processed_table/{}_rc.txt'.format(sample), sep ='\t').set_index('var')
    return(df_rc)

def import_lce(sample):
    df_lce = pd.read_csv('processed_table/{}_lce.txt'.format(sample), sep ='\t').set_index('var')
    return(df_lce)

def import_homo():
    df_homo =  pd.read_csv('processed_table/homo.txt', sep ='\t', index_col = 0)
    return(df_homo)

def combine_all():
    # wt = 'Variant_006995'
    # df_all.loc['Variant_004622']
    
    df_all = import_struct_with_classification()
    df_gce = import_gce()
    df_gce.columns = ['gce_dro', 'gce_mp']
    df_homo = import_homo()
        
    df_lce_dro = import_lce('dro')
    df_lce_mp = import_lce('mp')
    df_rc_dro = import_rc('dro')
    df_rc_mp = import_rc('mp')

    df_all = pd.merge(df_all, df_gce, left_index=True, right_index = True, how = 'inner')
    df_all = pd.merge(df_all, df_homo, left_index=True, right_index = True, how = 'outer')

    for i in range(-3, 4):
        df_all[f'rc{str(i)}_dro'] = df_rc_dro[str(i)]
    
    for i in range(-1, 3):
        df_all[f'lce{str(i)}_dro'] = df_lce_dro[str(i)]

    for i in range(-3, 4):
        df_all[f'rc{str(i)}_mp'] = df_rc_mp[str(i)]

    for i in range(-1, 3):
        df_all[f'lce{str(i)}_mp'] = df_lce_mp[str(i)]

    df_all.columns = ['main_loop', 'name_struct', 'seq', 'sg', 'motif', 'type_struct',
           'short_name', 'short_name_from_rand', 'pos1', 'pos2', 'pos3', 'pos4',
           'pos5', 'pos6', 'gce_dro', 'gce_mp', 'mcl_mp', 'homo_mp', 'mcl_dro', 'homo_dro',
           'rc-3_dro', 'rc-2_dro', 'rc-1_dro', 'rc0_dro', 'rc1_dro', 'rc2_dro',
           'rc3_dro', 'lce-1_dro', 'lce0_dro', 'lce1_dro', 'lce2_dro', 'rc-3_mp',
           'rc-2_mp', 'rc-1_mp', 'rc0_mp', 'rc1_mp', 'rc2_mp', 'rc3_mp',
           'lce-1_mp', 'lce0_mp', 'lce1_mp', 'lce2_mp']

    df_all.to_csv('processed_table/all.txt', sep ='\t')
    
    #df_all['struct'].value_counts()
    return(df_all)

def import_all():
    df_all = pd.read_csv('processed_table/all.txt', sep ='\t', index_col =0)
    df_all['short_name'] =df_all['short_name'].fillna('')
    df_all['short_name_from_rand'] =df_all['short_name_from_rand'].fillna('')
    return(df_all)


#%% (DONE) CALCULATE MOTIF COMBINATION

def make_combinations_for_motifs_each_structure(s):
    '''
    make the combination of base-pair or mismatches in each structures
    df_all_symm = import_df_all_symm('dro')
    df_all_symm = import_df_all_symm('14')

    '''
    df_all = import_all()
    
    # s = 'dro'
    df_all_symm = df_all[df_all['type_struct'] == 'symm'].copy() #mm: mismatch
    df_all_symm = df_all_symm[df_all_symm[f'rc0_{s}'] >= 0]
    df_all_symm['name_struct'].value_counts()
    df_all_symm = df_all_symm[df_all_symm.groupby(['main_loop'])['main_loop'].transform('count') >= 10]

    df_comb_all = pd.DataFrame()

    for structure in df_all_symm['short_name_from_rand'].unique():
        print(structure)
        # structure = '02-SSSSS SSSSS-02; '
        # structure = '01-SS SS-01; 04-S S-04; 06-S S-06; '
        df_comb = pd.DataFrame()
        
        df_each_struct = df_all_symm[df_all_symm['short_name_from_rand'] == structure].copy()
        
        dict_set_bp_each_position = {}
        for pos in [f'pos{i}' for i in range(1,7)]:
            dict_set_bp_each_position[pos] = tuple(df_each_struct[pos].unique().tolist())
        
        for number_nt in [1,2,3,4,5]:        
            for list_pos in itertools.combinations([f'pos{i}' for i in range(1,7)], number_nt):
                arrays = [dict_set_bp_each_position[pos] for pos in list_pos]
                # print(list(list_pos), arrays)
                
                for list_bp in list(itertools.product(*arrays)):
                    # print(list(list_bp))
                
                    df_using = df_each_struct.copy()
                    name_index = ''
    
                    for pos, bp in zip(list(list_pos), list(list_bp)):
                        df_using = df_using[df_using[pos] == bp]
                        name_index = name_index + f'{pos} {bp}; ' 
                    df_using = df_using[[f'homo_{s}',
                                         f'rc-3_{s}', f'rc-2_{s}', f'rc-1_{s}', 
                                         f'rc0_{s}', f'rc1_{s}', f'rc2_{s}', 
                                         f'rc3_{s}',  f'lce0_{s}',
                                         ]].dropna().copy()
                    
                    if len(df_using) > 0:
                        # count += 1
                        df_comb[name_index] = [len(df_using), df_using[ f'homo_{s}'].mean(),
                                           df_using[f'rc-3_{s}'].mean(), 
                                           df_using[f'rc-2_{s}'].mean(), df_using[f'rc-1_{s}'].mean(),
                                           df_using[f'rc0_{s}'].mean(), df_using[f'rc1_{s}'].mean(), 
                                           df_using[f'rc2_{s}'].mean(), df_using[f'rc3_{s}'].mean(),
                                           df_using[f'lce0_{s}'].mean(),
                                           structure]
                        
        df_comb = df_comb.T.reset_index()
        df_comb.columns = ['comb', f'number_variant_{s}', f'homo_{s}',
                           f'rc-3_{s}', f'rc-2_{s}', f'rc-1_{s}', 
                                         f'rc0_{s}', f'rc1_{s}', f'rc2_{s}', 
                                         f'rc3_{s}',  f'lce0_{s}', 'structure']
    
        df_comb_all = pd.concat([df_comb_all, df_comb], axis = 0)
        
    return(df_comb_all.set_index(['comb', 'structure']))
        

def merge_comb():
    df_dro = make_combinations_for_motifs_each_structure('dro') 
    df_mp = make_combinations_for_motifs_each_structure('mp') 
    
    df_dro.reset_index().to_csv('dro_comb.txt', sep = '\t', index = False)
    df_mp.reset_index().to_csv('mp_comb.txt', sep = '\t', index = False)
    
    df_comb_all  = pd.concat([df_dro, df_mp], axis =1,).reset_index()
    
    for pos in [f'pos{i}' for i in range(1,7)]:      
        df_comb_all[pos] = df_comb_all['comb'].map(lambda x: '' if pos not in x else 
                                                [i for i in x.split('; ') if pos in i][0].replace(f'{pos} ', ''))
    
    # df_comb_all = df_comb_all.reset_index()
    # df_comb_all = df_comb_all.drop('index', axis = 1)
    df_comb_all['group'] = df_comb_all['comb'].map(lambda x: '{} motif'.format(x.count(';')))
    
    df_comb_all.to_csv('motif_combination_rc_all_structures.txt', sep ='\t', index = False)

    df_comb_6bp = df_comb_all[df_comb_all['structure'] == '']
    df_comb_6bp.to_csv('./processed_table/motif_combination_rc_6bp.txt', sep ='\t', index = False)


def import_motif_combination():
    df_comb = pd.read_csv('./processed_table/motif_combination_rc_6bp.txt', sep ='\t', low_memory=False)
    for col in ['structure', 'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6']:
        df_comb[col] = df_comb[col].fillna('')
    
    # df_comb = pd.read_csv(f'./processed_table/motif_combination_rc_all_structures.txt', sep ='\t',).fillna('')
    return(df_comb)

# def import_best_motif():
#     path ='C:/DUC_research/Lab project/Non canonical mechanism/lib3_20211123/2.analyse/processed_table'
#     df_motif = pd.read_csv(f'{path}/top47_motifs.txt', sep ='\t').fillna('').set_index('motif_name')
#     return(df_motif)

#%% (DONE) MULTIPLE REGRESSION MODEL

def fit_all_36_2bp_combinations():
    df_all = import_all()
    df_all_6bp = df_all[df_all['name_struct'] == 'struct_0001']
    df_using = df_all_6bp[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'rc0_dro']].dropna().copy()

    X = df_using[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6']].copy()
    for list_pos in itertools.combinations([f'pos{i}' for i in range(1,7)], 2):
        X[' '.join(list_pos)] = X[list_pos[0]] +' '+ X[list_pos[1]]

    X = pd.get_dummies(data=X)
    Y = df_using['rc0_dro']

    model = sm.OLS(Y, X).fit()
    model.rsquared_adj
    print(model.summary())

    df_fit = pd.read_html(model.summary().tables[1].as_html(), header=0, index_col=0)[0]
    df_fit.columns = ['coef', 'std err', 't', 'p-val', 'ci_0.025', 'ci_0.975']

    df_fit.to_csv('./processed_table/fit all combinations rc0_dro.txt', sep = '\t', index = True, header = True)

    return(df_fit)

def import_fit_combinations():
    df_fit = pd.read_csv('./processed_table/fit all combinations rc0_dro.txt',sep = '\t', index_col = 0)
    return(df_fit)

#%% (DONE) CALCULATE SHIFTING SCORE (NO MISMATCH ONLY)

'''
Pos 0: C-G
Pos 7: T-A

given a combination: 4 nt
For example:
    pos1 A-T; pos2 A-T; pos3 A-T; pos4 T-A
    
    move motif: 
        --> pos1 A-T; pos2 A-T; pos3 A-T; pos4 T-A
'''

def return_comb_moved(comb):
    
    # comb = 'pos1 C-G; pos2 A-T; pos4 G-C; pos5 C-A; '
    # start_postion = [int(s) for s in re.findall(r'\b\d+\b', comb.replace('pos', ''))][0]
    # end_postion = [int(s) for s in re.findall(r'\b\d+\b', comb.replace('pos', ''))][-1]
    
    list_bp = [i.split(' ')[-1] for i in comb.split('; ')][:-1]
    list_pos = [int(s) for s in re.findall(r'\b\d+\b', comb.replace('pos', ''))]

    list_comb_shift = []    
    for i in range(-4, 4):
        if (i+ min(list_pos) >= 0) and (i+max(list_pos) <= 7):
            
            if ([i + x for x in list_pos][0] == 0) and (list_bp[0] != 'C-G'):
                continue
            elif ([i + x for x in list_pos][-1] == 7) and (list_bp[-1] != 'T-A'):
                continue
            else:
                comb_shift = '' 
                for pos, motif in zip([i + x for x in list_pos], list_bp):
                    comb_shift = comb_shift + f'pos{pos} {motif}; '    
                list_comb_shift.append(comb_shift.replace('pos0 C-G; ', '').replace('pos7 T-A; ', ''))
    return(list_comb_shift)

def calculate_shift_score(sample):
    '''
    calculate shited score
    '''
    # sample = 'dro'
    sample = 'mp'
    
    df_comb = import_motif_combination()
    df_comb['structure'].value_counts()
    df_comb['group'].value_counts()
    df_comb = df_comb[df_comb['structure'] == ''].set_index('comb') # select only match structure

    df_comb.columns
    
    
    df_shift_score = df_comb[[f'number_variant_{sample}',
                                     f'homo_{sample}', f'rc-3_{sample}', f'rc-2_{sample}',
                                     f'rc-1_{sample}', f'rc0_{sample}', f'rc1_{sample}',
                                     f'rc2_{sample}', f'rc3_{sample}', f'lce0_{sample}', 
                                     'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'group']]
    
    df_shift_score = df_shift_score[(df_shift_score[f'number_variant_{sample}'] >= 5) & 
                             (df_shift_score['group'].isin(['1 motif', '2 motif', '3 motif', '4 motif', '5 motif']))].copy()

    shift_score_list = []
    shift_comb_count = []
    for comb in df_shift_score.index:
        # comb = 'pos1 A-T; pos2 A-T; pos3 A-T; pos4 T-A; pos5 T-G; '
        list_comb_shift = return_comb_moved(comb)
        
        list_csites = [i - list_comb_shift.index(comb) for i in range(0, len(list_comb_shift))]
    
        list_rc = []
        list_rc_vs_homo = []
        for comb_shift, csites  in zip(list_comb_shift, list_csites):
            if (comb_shift in df_shift_score.index) and (csites in range(-3, 4)):
                rc_x = df_shift_score[f'rc{csites}_{sample}' ][comb_shift]
                homo_x = df_shift_score[f'homo_{sample}'][comb_shift]
                
                list_rc.append(rc_x)
                list_rc_vs_homo.append(rc_x/homo_x)
                # list_rc_vs_homo.append(df_comb['rc2_cl' + str(csites)][comb_shift] / )
        
        # shift_score_list.append(np.mean(list_rc_vs_homo)*((len(list_rc)-1)**0.5) / len(list_rc))
        shift_score_list.append(np.mean(list_rc_vs_homo) * np.mean(list_rc) *((len(list_rc))**0.5))
        shift_comb_count.append(len(list_rc))
    
    df_shift_score['shift_score'] = shift_score_list
    # df_using['ratio_score'] = ratio_score_list
    df_shift_score['shift_comb_count'] = shift_comb_count

    df_shift_score.to_csv(f'processed_table/rc_{sample}_shifting_score_nomismatch.txt', sep ='\t',)
    return()

def import_shift_score(sample):
    df_shift_score = pd.read_csv(f'./processed_table/rc_{sample}_shifting_score_nomismatch.txt', sep = '\t',
                                 index_col = 0, low_memory=False)
    for col in [ 'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6']:
        df_shift_score[col] = df_shift_score[col].fillna('')
    return(df_shift_score)


#%%

def make_dres():
    # using enrichment and shifting score
    df_using = pd.read_csv('./processed_table/topmotif_enrichment.txt', 
                                 sep = '\t', index_col = 0).fillna('')
    
    df_using = df_using[df_using['count_top'] >= 3]

    df_shift_dro = import_shift_score('dro')
    df_using['shift_score_dro'] = df_shift_dro['shift_score']
    df_using['rc0_dro'] = df_shift_dro['rc0_dro']
    df_using['shift_comb_count_dro'] = df_shift_dro['shift_comb_count']
    
    df_shift_mp = import_shift_score('mp')
    df_using['shift_score_mp'] = df_shift_mp['shift_score']
    df_using['rc0_mp'] = df_shift_mp['rc0_mp']
    df_using['shift_comb_count_mp'] = df_shift_mp['shift_comb_count']
    
    
    df_using = df_using[df_using['shift_score_dro'] >= 0.6]
    df_using = df_using[df_using['shift_comb_count_dro'] > 1]
    df_using = df_using[df_using['enrichment'] >= 2]
    

    df_dres = df_using.copy()
    df_dres['motif'] = df_dres.apply(lambda x: make_motif(x['pos1'], x['pos2'],x['pos3'], 
                                                          x['pos4'], x['pos5'], x['pos6']), axis = 1)
    
    df_dres['dres_score'] = (df_dres['shift_score_dro'] + df_dres['shift_score_mp'] ) /2
    
    df_dres.columns
    df_dres = df_dres.sort_values('dres_score', ascending = False)

    df_dres['motif_name'] = ['motif ' + str(i) for i in range(1, len(df_dres)+1)]
    df_dres['group'].value_counts()
    df_dres.to_csv('./processed_table/list_dres.txt', sep = '\t', index =True)


def import_dres():
    df_dres = pd.read_csv('./processed_table/list_dres.txt',sep ='\t', index_col = 0) 
    for i in range(1,7):
        df_dres[f'pos{i}'] = df_dres[f'pos{i}'].fillna('')
    return(df_dres)

