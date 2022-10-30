# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:51:05 2021

@author: Trung Duc Nguyen
"""

#%% IMPORT LIBRARY

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import using_now as us 
import common_tools as ct
import itertools
import statsmodels.api as sm

#%% IMPORT DATA

df_rc_mp = us.import_rc('mp')
df_rc_dro = us.import_rc('dro')

df_lce_mp = us.import_lce('mp')
df_lce_dro = us.import_lce('dro')

df_struct = us.import_struct_with_classification()
df_all = us.import_all()
df_all_6bp = df_all[df_all['name_struct'] == 'struct_0001']
df_comb_6bp = us.import_motif_combination()

df_shift_dro = us.import_shift_score('dro')
df_shift_mp = us.import_shift_score('mp')
df_dres = us.import_dres()
cutoff_shift = 0.6
cutoff_enrichment = 2

#%% (TOOL-DONE) HEATMAP CLEAVAGE ACCURACY ALL VARIANTS 

def draw_heatmap_rc_selected_motif():
    
    df_using = pd.read_csv('./processed table/topmotif_enrichment.txt', sep = '\t').fillna('')
    df_using = pd.merge(df_using, df_shift_dro[['shift_score']], left_on = 'comb', right_index= True)
    df_using['motif'] = df_using.apply(lambda x: us.make_motif(x['pos1'], x['pos2'],x['pos3'], 
                                                          x['pos4'], x['pos5'], x['pos6']), axis = 1)
    
    df_using= df_using[df_using['group'] == '4 motif']
    df_using = df_using[df_using['motif'] == 'CANGTN']

    # motif rank 1
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [1, 2, 3, 4, 5],
                                            [['C-G'], ['A-T'], ['T-A'], ['G-C'], ['A-T']],
                                            'motif #1 rc0_dro ', draw = True)

    # motif rank 47
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [ 2, 4, 5, 6],
                                            [ ['A-T'], ['G-C'], ['A-T'], ['G-T']],
                                            'motif #47 rc0_dro ', draw = True)

    # motif rank 173
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [ 2, 4, 6],
                                            [ ['A-T'], ['G-C'], ['G-T']],
                                            'motif #173 rc0_dro ', draw = True)

    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [ 1, 3, 5],
                                            [ ['A-T'], ['G-C'], ['G-T']],
                                            'motif #173_-1 rc0_dro ', draw = True)

    # nonDRES
    '''
        df_dres = us.import_dres()
    
        list_3bp_dres = df_dres[df_dres['group'] == '3 motif'].index.tolist()
        list_4bp_dres = df_dres[df_dres['group'] == '4 motif'].index.tolist()
        list_5bp_dres = df_dres[df_dres['group'] == '5 motif'].index.tolist()
    
        df_all_dres = pd.DataFrame()
    
        for motif in list_3bp_dres + list_4bp_dres +  list_5bp_dres:
            df =  us.extract_variant_given_motif(df_all_6bp, motif)
            df_all_dres = pd.concat([df, df_all_dres], axis = 0)
    
        df_nondres = df_rc_dro[~df_rc_dro.index.isin(df_all_dres.index)]
        df_nondres = df_nondres[df_nondres.index.isin(df_all_6bp.index)]
    
        us.draw_heatmap_cleavage_accuracy(df_nondres.copy(), "heatmap rc_dro nondres")
    
    '''

    # motif B
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [1, 2, 3, 4,],
                                            [['A-T'], ['C-G'], ['G-C'], ['T-A']],  'motif B1 rc0_dro ')
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [2, 3, 4, 5],
                                            [['A-T'], ['C-G'], ['G-C'], ['T-A']],  'motif B2 rc0_dro ')
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [ 3, 4, 5, 6],
                                            [['A-T'], ['C-G'], ['G-C'], ['T-A']],  'motif B3 rc0_dro ')
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_dro, [4, 5, 6],
                                            [['A-T'], ['C-G'], ['G-C']],  'motif B4 rc0_dro ')
    
    # test
    '''
    pos2 G-C; pos3 G-C; pos4 G-C; pos5 T-A; pos6 G-T; 
    '''
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_mp, [2, 3, 4,5,6],
                                            [['G-C'], ['G-C'], ['G-C'], ['T-A'], ['G-T']], 
                                            'motif T2 rc0_mp ', draw = True)
    us.draw_heatmap_rc_pattern_given_pos_bp(df_all_6bp, df_rc_mp, [1,2, 3, 4,5],
                                            [['G-C'], ['G-C'], ['G-C'], ['T-A'], ['G-T']], 
                                            'motif T1 rc0_mp ', draw = True)

    
    # mp
    us.draw_heatmap_cleavage_accuracy(df_rc_mp.copy(), "heatmap rc all mp")
    
    us.draw_heatmap_cleavage_accuracy(df_rc_mp[df_rc_mp.index.isin(df_all[df_all['type_struct'] == 'symm'].index)],
                                      "heatmap rc symm struct mp")
    
    us.draw_heatmap_cleavage_accuracy(df_rc_mp[df_rc_mp.index.isin(df_all[df_all['name_struct'] == 'struct_0001'].index)],
                                      "heatmap rc 6-bp mp")
    
    # dro
    us.draw_heatmap_cleavage_accuracy(df_rc_dro.copy(), "heatmap rc all dro")
    us.draw_heatmap_cleavage_accuracy(df_rc_dro[df_rc_dro.index.isin(df_all[df_all['type_struct'] == 'symm'].index)],
                                      "heatmap rc symm struct dro")
    
    us.draw_heatmap_cleavage_accuracy(df_rc_dro[df_rc_dro.index.isin(df_all[df_all['name_struct'] == 'struct_0001'].index)],
                                      "heatmap rc 6-bp dro")
    
    
#%% (DONE) FIGURE 3: NUMBER OF VARIANT RECOVERED

def number_variant_recovered():
    '''
    number of variant >= 20
    '''
    
    df_using = us.import_sg()
    df_set2 = us.import_ctrl_count('set2')
    df_set1 = us.import_ctrl_count('set1')
    df_ctrl = pd.concat([df_set2, df_set1]).groupby(level=0).max()
    
    df_using['count'] = df_ctrl['count']

    df_using = df_using.dropna()
    df_using = df_using.groupby('sg')['motif'].count().to_frame()
    
    df_using['motif'].sum()

    fig, ax = plt.subplots(1, 1, figsize=(1,3))
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})

    sns.stripplot(y = df_using['motif'], orient='v', color = 'white', linewidth = 0.5,  jitter=0.25,
                  ec ='k', fc ='whitesmoke', alpha = 0.6,  s = 7)
    plt.ylim((3000, 4200))
    plt.yticks([3000, 3500,  4096], fontsize = 0, color ='white',  )
    sns.despine()
    plt.ylabel('',fontsize =0, color = 'white')

    plt.xlabel('',fontsize =0, color = 'white')

    ct.decorate_boxplot(ax)
    plt.xticks([], fontsize = 0, color ='white',)
    # plt.yticks(fontsize = 0, color ='white',)
    plt.savefig('stripplot variant recovered', dpi=300, bbox_inches="tight") 

#%% (DONE) FIGURE 3: STACKED BARPLOT NUMBER OF STRUCTURES RECOVERED

def figure_number_structure_recovered():
    df_using = us.import_all()
    # df_using = df_using[df_using.groupby(['main_loop'])['main_loop'].transform('count') >= 20]
    df_using = df_using.groupby(['main_loop', 'type_struct', 'short_name_from_rand'])['short_name'].count().reset_index()

    df_using['type_struct'].value_counts()
    df_using.groupby('type_struct')['short_name'].sum() / 253679

    # draw bar using illustrator
    '''

    type_struct
    asymm      0.016655
    symm       0.735272
    unknown    0.248073
    
    '''
#%% (DONE) FIGURE 3: BARPLOT NUMBER OF BASE PAIRS RECOVERED

def draw_barplot_number_variants_with_different_bp():
    
    df_struct_sym = df_struct[df_struct['type_struct'] == 'symm'].copy()
    df_struct_sym['number_mismatch'] = df_struct_sym['short_name'].map(lambda x: x.count('S')/2)
    list_expected_bp_variants = df_struct_sym['number_mismatch'].value_counts().to_frame().sort_index()['number_mismatch'].tolist()
    
    sum(list_expected_bp_variants)

    # calculate actual variant recoverved
    df_using = df_all[df_all.index.isin(df_struct_sym.index)].copy()
    list_actual_bp_variants = []
    
    list_basepair = ['A-T', 'T-A', 'C-G', 'G-C', 'T-G', 'G-T']
    df_using['pos1'] = df_using['motif'].map(lambda x: 0 if x[0]+ '-' +x[12] in list_basepair else 1)
    df_using['pos2'] = df_using['motif'].map(lambda x: 0 if x[1]+ '-' +x[11] in list_basepair else 1)
    df_using['pos3'] = df_using['motif'].map(lambda x: 0 if x[2]+ '-' +x[10] in list_basepair else 1)
    df_using['pos4'] = df_using['motif'].map(lambda x: 0 if x[3]+ '-' +x[9] in list_basepair else 1)
    df_using['pos5'] = df_using['motif'].map(lambda x: 0 if x[4]+ '-' +x[8] in list_basepair else 1)
    df_using['pos6'] = df_using['motif'].map(lambda x: 0 if x[5]+ '-' +x[7] in list_basepair else 1)
    df_using['number_mismatch'] = df_using[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6']].sum(axis = 1)
    
    list_actual_bp_variants = df_using['number_mismatch'].value_counts().to_frame().sort_index()['number_mismatch'].tolist()
    
    # draw 
    fig, ax = plt.subplots(1, 1, figsize=(4,4))
    
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})

    ax.bar(range(6), list_expected_bp_variants[::-1], align='center', 
            color = 'lightgrey',
            linewidth = 2, width = 0.8)
    
    ax.bar(range(6), list_actual_bp_variants[::-1], align='center', 
            color ='#60d394', linewidth = 2, width = 0.8)
    
    ct.decorate_boxplot(ax)
    sns.despine()
    
    plt.ylim([0,100000])
    plt.yticks([0, 20000, 40000, 60000, 80000, 100000],  color ='white', size =0)
    
    plt.xticks([0,1,2,3,4,5],  color ='white', size =0)
    
    plt.savefig('stackedbarplot number of variants recovered with different number of base pairs',
                bbox_inches="tight", dpi = 300)

#%% (DONE) FIGURE 3: BARCODE DISTRIBUTION

def histogram_barcode_distribution():
    
    df_set2 = us.import_ctrl_count('set2')
    df_set1 = us.import_ctrl_count('set1')
    df_ctrl = pd.concat([df_set2, df_set1]).groupby(level=0).max()
    
    # df_ctrl = df_ctrl[df_ctrl['count'] >= 20]
    
    fig, ax = plt.subplots(1, 1, figsize=(4,4))
    
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})
    ct.decorate_boxplot(ax)


    ax.hist( np.log10(df_ctrl['barcode_number']), histtype='step', 
             bins = 100, linewidth =2, color ='grey')

    print(df_ctrl.mean()) # 197.665822
    print(df_ctrl.median()) # 127.0

    plt.ylim([0,10000])
    plt.xlim([0,4])
    
    plt.xticks([0, 1, 2,3, 4], color='white', size =0)
    plt.yticks([0,2000,  4000, 6000, 8000, 10000], color='white', size =0)
    
    sns.despine()
    plt.savefig('histogram barcode number distribution',  bbox_inches="tight", dpi =300)

#%% (DONE) FIGURE 3: DRAW LINE PLOT MAJOR CLEAVAGE SITE

def draw_line_plot_major_cleavage_site():
    
    df_draw = pd.DataFrame()
    
    df_using = df_all.copy()
    # df_using = df_all[df_all['name_struct'] == 'struct_0001'].copy() # for 6bp
    df_using = df_all[df_all['type_struct'] == 'symm'].copy() # for symm struct
    
    df_draw['dro'] = df_using['mcl_dro'].dropna().astype(int).value_counts()
    df_draw['mp'] = df_using['mcl_mp'].dropna().astype(int).value_counts()
    
    df_draw = df_draw.sort_index()
    df_draw = df_draw / df_draw.sum()*100
    
    # draw
    fig, ax = plt.subplots(1, 1, figsize=(5,2))
    
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})

    ct.decorate_boxplot(ax)

    ax.plot(range(1,12), df_draw['mp'], marker = 'o', ms=10, color = '#60d394', lw =2)
    ax.plot(range(1,12), df_draw['dro'], marker = 'o', ms=10, color = '#3da4dc', lw =2)
    
    sns.despine()
    plt.xticks(range(1,12), fontsize =0, color = 'white')
    plt.ylim(-15,80)
    plt.yticks([ 0, 40, 80],  fontsize =0, color = 'white')
    
    plt.ylabel('',fontsize =0, color = 'white')
    
    plt.xlabel('',fontsize =0, color = 'white')
    
    # plt.savefig('lineplot fraction of variants with different major cleavage site all variants',
    plt.savefig('lineplot fraction of variants with different major cleavage site symmetric variants',
    # plt.savefig('lineplot fraction of variants with different major cleavage site 6-bp variants',
                bbox_inches="tight",dpi = 300)

#%% (DONE) FIGURE 3: DRAW GCE DRO MP

def boxplot_compare_gce_dro_mp():
    df_using = df_all.copy()
    
    df_using = df_using[df_using['name_struct'] == 'struct_0001']
    # df_using = df_using[df_using['type_struct'] == 'symm']
    df_using.columns
    df_draw = df_using[[ 'gce_dro', 'gce_mp']].melt()

    # draw
    fig, ax = plt.subplots(1, 1, figsize=(2,4))
    
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})

    boxprops = {'edgecolor': 'k', 'linewidth': 1.5}
    lineprops = {'color': 'k', 'linewidth': 1.5}
    boxplot_kwargs = {'boxprops': boxprops,  'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,}

    sns.boxplot(data = df_draw, x = 'variable', y = 'value', palette = ['#3da4dc', '#60d394'],
                fliersize=1, showfliers=True, **boxplot_kwargs, saturation=1)

    ct.decorate_boxplot(ax)
    ct.adjust_box_widths(fig, 0.8)

    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    sns.despine()
    plt.ylim(-12, 4)
    
    plt.xticks( fontsize =0, color = 'white')
    plt.xlabel('',fontsize =0, color = 'white')

    plt.yticks( fontsize =0, color = 'white') 
    plt.ylabel('',fontsize =0, color = 'white')
    # plt.savefig( 'boxplot gce all variants dro and mp', bbox_inches= 'tight', dpi= 300)
    # plt.savefig( 'boxplot gce all symmetric variants dro and mp', bbox_inches= 'tight', dpi= 300)
    plt.savefig( 'boxplot gce 6-bp variants dro and mp', bbox_inches= 'tight', dpi= 300)
            
    # statistic
    ct.statistic_mannwhitneyu(df_using['gce_dro'].dropna(), df_using['gce_mp'].dropna()) # '0.0'

 

#%% (DONE) FIGURE 3: CLEAVAGE SITE CHANGE DROSHA AND MP

def cleavage_site_change_drosha_mp():

    df_using = df_all.copy()
    df_using = df_using[df_using['type_struct']=='symm']    
    df_using = df_using[['mcl_dro', 'mcl_mp', 'rc0_dro','rc0_mp',  'homo_dro', 'homo_mp']]
    df_using = df_using.dropna()
    
    df_using['mcl_consistence'] = df_using.apply(lambda x: 'consistent' if x['mcl_dro'] == x['mcl_mp'] else
                                                 'inconsistent',axis = 1)
    df_using['mcl_consistence'].value_counts()

    def example_for_mcl():
        df_using2 = df_using[df_using['mcl_consistence'] == 'consistent'].copy()
        df_using2 = df_using2[df_using2['rc0_mp'] >= 0.9]
        df_using2 = df_using2[df_using2['rc0_dro'] >= 0.4]
        
        # Variant_131444
        df_rc_dro.loc['Variant_131444']
        df_rc_mp.loc['Variant_131444']
        
        variant = 'Variant_131444'
        sample = 'mp'
        sample = 'dro'
        
        fig, ax = plt.subplots(1, 1, figsize=(5,3))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': False,})
        
        if sample == 'dro':
            plt.bar(range(0,7), df_rc_dro[[str(i)  for i in range(-3,4)]].loc[variant], 
            # plt.bar(range(0,7), df_rc_mp.loc['Variant_077194'], 
                    color = ['lightgrey']*3 + ['#3da4dc'] +['lightgrey']*3)
        else:        
            plt.bar(range(0,7), df_rc_mp[[str(i) for i in range(-3,4)]].loc[variant], 
                    color = ['lightgrey']*3 + ['#3da4dc'] +['lightgrey']*3)
        plt.ylim(0,1)
        sns.despine()
        ct.decorate_boxplot(ax)
        plt.ylabel('',fontsize =0, color = 'white')
    
        plt.xlabel('',fontsize =0, color = 'white')
        plt.xticks(range(0,7), fontsize = 0, color ='white',)
        plt.yticks([0,1], fontsize = 0, color ='white',)
        
        plt.savefig(f'barplot {sample} {variant} example for mcl', bbox_inches ='tight', dpi = 300)
    
    # pie plot consistent and inconsistent mCL
    
    def pieplot_consistent_inconsistent_mCL(): # DONE
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
        
        # 65538 / 114699 = 57.14%
        # 49161 / 114699 = 42.86%
        plt.pie([65538, 49161],startangle=0, colors=['#3da4dc', 'lightgrey', ], 
                explode=(0, 0.1 ), wedgeprops={"edgecolor":"white"})
        
        plt.savefig('pieplot consistent and inconsistent mcl between drosha and mp.png', bbox_inches="tight", dpi = 300)

    def densityplot_homo_consistent_variants(): # DONE
    
        type_ = 'consistent'
        type_ = 'inconsistent'
        # dgcr8 further enhance the cleavage of Drosha
        # df_draw = df_using[df_using['mcl_consistence'] == 'consistent'].copy()
        df_draw = df_using[df_using['mcl_consistence'] == type_].copy()
        
        print(len(df_draw[df_draw['homo_dro'] > df_draw['homo_mp']]))
        print(len(df_draw[df_draw['homo_dro'] <= df_draw['homo_mp']]))
        
        fig, ax = plt.subplots(1, 1, figsize=(4,4))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
    
        if type_ == 'consistent':
            sns.kdeplot(data= df_draw, x= 'homo_dro', y = 'homo_mp', levels=5, alpha = 0.7,
                        fill=True, color = '#3da4dc')
        
        else :
        
            sns.kdeplot(data= df_draw, x= 'homo_dro', y = 'homo_mp', levels=5, alpha = 0.7,
                        fill=True, color = 'lightgrey')
        
        plt.plot([0,1], [0,1], lw = 1, ls = '--', c = 'k')
        sns.despine()
        
        plt.ylabel('',fontsize =0, color = 'white')
    
        plt.xlabel('',fontsize =0, color = 'white')
    
        ct.decorate_boxplot(ax)
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 0, color ='white',)
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 0, color ='white',)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.savefig(f'density plot homo of {len(df_draw)} {type_} variants', bbox_inches='tight', dpi = 300)

    def lineplot_number_variants_with_different_mcl_consistent ():
        # dgcr8 made drosha cleave at cl0
        # line plot number of variants and mCL
        
        # df_using2 = df_using[df_using['mcl_consistence'] == 'inconsistent'].copy()
        df_using2 = df_using[df_using['mcl_consistence'] == 'consistent'].copy()
    
        
    
        df_draw = pd.DataFrame()
        df_draw['dro'] = df_using2['mcl_dro'].value_counts()        
        df_draw['mp'] = df_using2['mcl_mp'].value_counts()        
        
        df_draw = df_draw.reindex([i for i in range(-5, 6)])
    
        # draw
        fig, ax = plt.subplots(1, 1, figsize=(5,3))
        
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
    
        ct.decorate_boxplot(ax)
    
        # inconsistent
        # ax.plot(range(1,12), df_draw['mp'], marker = 'o', ms=10, color = '#60d394', lw =2)
        # ax.plot(range(1,12), df_draw['dro'], marker = 'o', ms=10, color = '#3da4dc', lw =2)
        # plt.ylim(-5000,40000)
        # plt.yticks([ 0, 10000, 20000, 30000, 40000],  fontsize =0, color = 'white')
        
        # consistent
        ax.plot(range(1,12), df_draw['mp'], marker = 'o', ms=10, color = '#60d394', lw =2, alpha = 0.5)
        ax.plot(range(1,12), df_draw['mp'], marker = 'o', ms=10, color = '#3da4dc', lw =2, alpha = 0.5)
        plt.ylim(-5000,50000)
        plt.yticks([ 0, 10000, 20000, 30000, 40000, 50000],  fontsize =0, color = 'white')
        
        sns.despine()
        plt.xticks(range(1,12), fontsize =0, color = 'white')
        plt.ylabel('',fontsize =0, color = 'white')
        
        plt.xlabel('',fontsize =0, color = 'white')
        
        plt.savefig('line plot number of 65538 consistent variants with different major cleavage site',
        # plt.savefig('line plot number of 49161 inconsistent variants with different major cleavage site',
                    bbox_inches="tight",dpi = 300)

        def lineplot_average_rc_inconsistent_variant():
            
            group = 'consistent'
            # group = 'inconsistent'
            value = 'rc'
            # value = 'lce' # lce not good to compare 2 samples because its related to normal 
            
            df_using2 = df_using[df_using['mcl_consistence'] == group].copy()
            
            if value == 'rc':
                df1 = df_rc_dro[df_rc_dro.index.isin(df_using2.index)].copy()
                df2 = df_rc_mp[df_rc_mp.index.isin(df_using2.index)].copy()
            else:
                df1 = df_lce_dro[df_lce_dro.index.isin(df_using2.index)].copy()
                df2 = df_lce_mp[df_lce_mp.index.isin(df_using2.index)].copy()
                
            df1['sample']= 'dro'
            df2['sample']= 'mp'

            df1 = df1.melt(id_vars ='sample')
            df2 = df2.melt(id_vars ='sample')

            df = pd.concat([df1, df2], axis = 0)

            fig, ax = plt.subplots(1, 1, figsize=(5,3))
            
            sns.set_style('ticks', {'axes.edgecolor': 'black',
                                    'grid.linestyle': 'dashed',
                                    'axes.grid': True,})
            
            ct.decorate_boxplot(ax)

            sns.lineplot(data = df, hue = 'sample', x = 'variable', y = 'value', 
                         marker = 'o', ms=12, err_style="band", palette =['#3da4dc', '#60d394', ] )
            # for rc
            if value == 'rc':
                plt.ylim(0,0.5)
                plt.yticks([ 0, 0.1, 0.2, 0.3, 0.4, 0.5],  fontsize =0, color = 'white')
            # else:
                # plt.ylim(-5,0)
                # plt.yticks([ -5, -4, -3, -2, -1, 0],  fontsize =0, color = 'white')
            # for lce
            # plt.ylim(0,0.5)
            # plt.yticks([ 0, 0.1, 0.2, 0.3, 0.4, 0.5],  fontsize =0, color = 'white')

            sns.despine()
            plt.xticks(range(0,11), fontsize =0, color = 'white')
            plt.ylabel('',fontsize =0, color = 'white')
            
            plt.xlabel('',fontsize =0, color = 'white')
            plt.legend([], frameon=False)
            
            plt.savefig(f'lineplot {value} {group} variants between dro and mp', bbox_inches = 'tight', dpi = 300)
    


#%% (DONE) FIGURE 4: ADJUSTED-R SQUARE DIFFERENT COMBINATION

def draw_bar_plot_compare_rsquare_different_combination():
    df_using = df_all_6bp.copy()
    
    # df_using = df_using[df_using['major_csite2'] == 0]
    # df_using = df_using[df_using['rc2_cl0'] >= 0]
    
    df_using = df_using[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'rc0_dro']].dropna()

    df_rsquared_adj = pd.DataFrame(columns = ['no_interact', 'interact'])
    
    Y = df_using['rc0_dro']
    
    for number_nt in [1,2,3,4,5,6]:        
        for list_pos in itertools.combinations([f'pos{i}' for i in range(1,7)], number_nt):
            X = df_using[list(list_pos)].copy()
            
            X1 = pd.get_dummies(data=X)
            if number_nt != 1:
                for list2_pos in itertools.combinations(list(list_pos), 2):
                    X[' '.join(list2_pos)] = X[list2_pos[0]] +' '+ X[list2_pos[1]]
            
            X2 = pd.get_dummies(data=X)
                
            model1 = sm.OLS(Y, X1).fit()
            model2 = sm.OLS(Y, X2).fit()
                
            df_rsquared_adj.loc[' + '.join(list(list_pos))] = [ model1.rsquared_adj, model2.rsquared_adj] 
            
            print(list_pos)

    df_rsquared_adj['number of position'] = df_rsquared_adj.index.map(lambda x: x.count('+') + 1)

    for pos in [f'pos{i}' for i in range(1,7)]:      
            df_rsquared_adj[pos] = df_rsquared_adj.index.map(lambda x: 
                                                    0 if pos not in x else 1)
        
    df1 = df_rsquared_adj[df_rsquared_adj['number of position'] == 1].copy().sort_values('interact')
    df2 = df_rsquared_adj[df_rsquared_adj['number of position'] == 2].copy().sort_values('interact')
    df3 = df_rsquared_adj[df_rsquared_adj['number of position'] == 3].copy().sort_values('interact')
    df4 = df_rsquared_adj[df_rsquared_adj['number of position'] == 4].copy().sort_values('interact')
    df5 = df_rsquared_adj[df_rsquared_adj['number of position'] == 5].copy().sort_values('interact')
    df6 = df_rsquared_adj[df_rsquared_adj['number of position'] == 6].copy().sort_values('interact')
    df0 = pd.DataFrame(data = {'no_interact': [0], 'interact': [0], 'number of position': [-1],
                               'pos1': [-1], 'pos2': [-1],
                               'pos3':[-1], 'pos4': [-1],
                               'pos5': [-1], 'pos6': [-1]} )
    
    df_draw = pd.concat([df1,df0, df2, df0, df3, df0, df4, df0, df5, df0, df6]).reset_index()
    
    
    df_draw['color'] = df_draw.apply(lambda x: '#60d394' if x['pos1'] + x['pos4'] == 2 else 'lightgrey', axis = 1)
    df_draw['color'] = df_draw.apply(lambda x: '#60d394' if x['index'] in ['pos1', 'pos4'] else x['color'], axis = 1 )
    
    # draw annotation
    fig, ax = plt.subplots(1, 1, figsize=(1,9))
    
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dotted',
                            'axes.grid': True,})
    
    for pos in range(1,7):
        for index in df_draw.index:
            if df_draw['pos' + str(pos)][index] == 1:
                ax.scatter(pos, index, s = 30,  color = 'lightgrey', 
                           edgecolors = 'k', linewidth=0.6,  alpha =0.75, zorder=3)
            elif df_draw['pos' + str(pos)][index] == 0:
                ax.scatter(pos, index, s= 30, c = 'white', edgecolors = 'grey', linewidth=0.2,)
                
    plt.yticks([], color = 'white', size =0) 
    plt.xticks([], color = 'white', size =0)
    plt.ylabel('',fontsize =0, color = 'white')
    plt.xlabel('',fontsize =0, color = 'white')
    
    sns.despine(left=True, right= True, bottom=True, top = True)
    plt.savefig('dot annotation of adjusted rsquare different position combinations.png', 
                bbox_inches="tight", dpi = 300)
            

    # draw bar plot
    fig, ax = plt.subplots(1, 1, figsize=(9, 1.5))
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dotted',
                            'axes.grid': True,})



    ax.bar(range(len(df_draw)), 
           df_draw['interact'], width=0.8,
           color =df_draw['color'].tolist() , 
           alpha = 1)
    
    # ax.bar(range(len(df_draw)), 
    #         df_draw['no_interact'], 
    #        # color ='#25a18e',
    #        color ='lightgrey',
    #        alpha = 1)
    
    
    ct.decorate_boxplot(ax)
    sns.despine()
    plt.xlim(-1, 68.5)
    plt.xticks([],  color ='white', size =0)
    plt.ylim(0, 0.8)
    plt.yticks([0,  0.2,  0.4,  0.6, 0.8],  color ='white', size =0)
    plt.savefig('barplot adjusted rsquare different position combinations fitting rc0_dro and position.png',
                bbox_inches="tight", dpi = 300)

#%% (DONE) FIGURE 4: CHORD DIADRAM - FIT ALL 36 + 2bp COMBINATIONS

# Calculate the xy coords for each point on the circle
def calculate_position_each_point(npoints = 36):
    s = 2 * np.pi / npoints
    list_x = []
    list_y = []
    
    for i in np.arange(npoints):
        angle = s * i
        x = npoints * np.cos(angle)
        y = npoints * np.sin(angle)
        list_x.append(x)
        list_y.append(y)
    return(list_x, list_y)

def draw_arc(pointA, pointB, color, linewidth,alpha, ax):
    # Calculate the centre of the Arc
    list_x, list_y = calculate_position_each_point(npoints = 36)
    
    x1 = list_x[pointA - 1]
    y1 = list_y[pointA - 1]
  
    x2 = list_x[pointB - 1]
    y2 = list_y[pointB - 1]

    bezier_path = np.arange(0, 1.01, 0.01)
    
    x = (1 - bezier_path)** 2 * x1 + 2 * (1 - bezier_path) * bezier_path * 0 + bezier_path** 2 * x2
    y = (1 - bezier_path)** 2 * y1 + 2 * (1 - bezier_path) * bezier_path * 0 + bezier_path** 2 * y2
    
    ax.plot(x, y, color, linewidth = linewidth, alpha = alpha, zorder = 1)


def draw_chord_diagram_coefficient():
    df_fit = us.import_fit_combinations()

    # process df1, df2
    df1 = df_fit[df_fit.index.str.count('pos') == 1].copy()
    df1 = df1.reindex(index = [x + '_' + y for x, y in itertools.product([f'pos{i}' for i in range(1,7)], 
                                                                             ['A-T', 'T-A', 'C-G', 'G-C', 'G-T', 'T-G'])])
    dict_position = {}
    for  bp, pos in zip(df1.index, range(1,37)):
        dict_position[bp] = pos
    
    df2 = df_fit[df_fit.index.str.count('pos') == 2].copy()
    df2['posA'] = df2.index.map(lambda x: x.split('_')[0].split(' ')[0] + '_' + x.split('_')[1].split(' ')[0])
    df2['posB'] = df2.index.map(lambda x: x.split('_')[0].split(' ')[1] + '_' + x.split('_')[1].split(' ')[1])
    
    df2['posA'] = df2['posA'].map(lambda x: dict_position[x])
    df2['posB'] = df2['posB'].map(lambda x: dict_position[x])
    df2['ab_coef'] = df2['coef'].map(lambda x: abs(x))
    
    # draw parameters
    # from  matplotlib.colors import LinearSegmentedColormap
    # cmap= LinearSegmentedColormap.from_list('rg',["darkred","red","lightcoral","white", "palegreen","green","darkgreen"][::-1],
    #                                         N=100) 
    
    vmax = 0.1
    vmin = -0.1
    cmap = 'RdYlGn_r'
    list_x, list_y = calculate_position_each_point(npoints = 36)
    fig, ax = plt.subplots(1, 1, figsize=(3,3))
        
    # draw connections
    for index in df2.sort_values('ab_coef').index:
        coef = df2['coef'][index] 
        if (coef >= 0.045) or (coef <= -0.045):
            a = df2['posA'][index]
            b = df2['posB'][index]
            color = ct.color_map_color(coef, cmap_name= cmap, vmin=vmin, vmax=vmax)
            draw_arc(a, b, color, abs(coef)*60, min(abs(coef)*6, 1),  ax) 
        
    # draw circles
    ax.scatter(list_x, list_y, marker='o', s=160, 
               c=[ ct.color_map_color(i, cmap_name=cmap, vmin=vmin, vmax=vmax) for i in df1['coef']], 
               alpha =1, zorder = 2, ec = 'k', lw = 1)
    
    sns.despine(top =True, bottom=True, left=True, right = True)
    plt.xticks([])
    plt.yticks([])
    
    plt.savefig('chord diagram fit all 36 bp and combination.png', bbox_inches='tight', dpi =300)


#%% (DONE) FIGURE 4: BOXPLOT INTERACTION EFFECT 

def draw_boxplot_interaction_effect(df_draw, save_fig):
    boxprops = {'edgecolor': 'k', 'linewidth': 1.5}
    lineprops = {'color': 'k', 'linewidth': 1.5}
    boxplot_kwargs = {'boxprops': boxprops,  'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,}

    
    fig, ax = plt.subplots(1, 1, figsize=(3,3))
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})
    ct.decorate_boxplot(ax)

    sns.boxplot(data = df_draw, x ='type', y = 'rc0_dro', order = ['', 'A', 'B', 'AB'],
                fliersize=0.5, showfliers=True, **boxplot_kwargs,
                palette=['white']*3 + ['#60d394'])

    ct.adjust_box_widths(fig, 0.8)

    plt.xticks( fontsize =0, color = 'white')
    plt.xlabel('',fontsize =0, color = 'white')

    plt.ylim(-0.05,1)
    plt.yticks([0,0.2, 0.4, 0.6, 0.8, 1.0], fontsize =0, color = 'white') 
    plt.ylabel('',fontsize =0, color = 'white')
    sns.despine()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    plt.savefig(save_fig,
                bbox_inches="tight",dpi = 300)

def run_draw_boxplot_interaction_effect():
    df_fit = us.import_fit_combinations()
    df_using = df_fit[df_fit.index.str.count('pos') == 2].copy()

    
    list_pair = df_using[abs(df_using['coef']) >= 0.05 ].index.tolist() + \
    df_using.sort_values('p-val', ascending = False).head(20).index.tolist()
    
    for pair in list_pair:
        if df_using['p-val'][pair]  > 0.7:
            effect = 'no effect'
        elif df_using['coef'][pair] > 0:
            effect = 'synergistic effect'
        elif df_using['coef'][pair] < 0:
            effect = 'antagonistic effect'
        
        
        '''
        pair = 'pos1 pos4_G-C A-T' # 0.0001 # no effect
        pair = 'pos2 pos4_A-T G-C' # 0.0681 # syner
        pair = 'pos3 pos4_G-C T-A' # -0.0572 # antago
        
        df_using['coef'][pair]
        df_using['p-val'][pair]
        
        '''
        
        
        posA = pair.split('_')[0].split(' ')[0]
        posB = pair.split('_')[0].split(' ')[1]
        bpA = pair.split('_')[1].split(' ')[0]
        bpB = pair.split('_')[1].split(' ')[1]
        
        df = df_all_6bp[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'rc0_dro']].dropna().copy()
        
        df['posA'] = df[posA].map(lambda x: 'A' if bpA in x else '' )
        df['posB'] = df[posB].map(lambda x: 'B' if bpB in x else '' )
        df['type'] = df['posA'] + df['posB']
        
        savefig = ' '.join([str(i) for i in df['type'].value_counts().tolist()])
        
        draw_boxplot_interaction_effect(df, 
                                        save_fig = f'Select boxplot {effect} {pair} {savefig} variants' )
    
        # calculate cohen d # non vs A, B vs AB
        print('A vs non', ct.calculate_cohend(df[df['type'] == 'A']['rc0_dro'].tolist(),
                                              df[df['type'] == '']['rc0_dro'].tolist()))
        
        print('AB vs B', ct.calculate_cohend(df[df['type'] == 'AB']['rc0_dro'].tolist(),
                                              df[df['type'] == 'B']['rc0_dro'].tolist()))   
        
        # calculate p-value
        print('A vs non', ct.statistic_mannwhitneyu(df[df['type'] == 'A']['rc0_dro'].tolist(),
                                              df[df['type'] == '']['rc0_dro'].tolist()))
        
        print('AB vs B', ct.statistic_mannwhitneyu(df[df['type'] == 'AB']['rc0_dro'].tolist(),
                                              df[df['type'] == 'B']['rc0_dro'].tolist()))   
        
        
#%% (DONE) FIGURE 4: CALCULATE ENRICHED SCORE AND MAKE DRES

def motif_generation_enriched_score_calculation():
    df_using = df_all_6bp.copy()
    df_using = df_using[['lce0_dro', 'rc0_dro', 'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6']]
    # df_using['oddratio_rc0_dro'] = df_using['rc0_dro'].map(lambda x: np.log2(x/(1-x))) 

    df_using = df_using.dropna()
    cutoff_rc0 = 0.5
    cutoff_lce = -4

    def draw_lce0_rc0_all_6bp_variants():
        # draw lce0 and rc0 all samples
        fig, ax = plt.subplots(1, 1, figsize=(4,4))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
        ct.decorate_boxplot(ax)
    
        sns.scatterplot(data = df_using,
                        x = 'rc0_dro', y = 'lce0_dro',
                        # x = 'oddratio_rc0_dro', y = 'lce0_dro',
                        linewidth = 0.1, alpha = 1,  s=1, color = 'lightgrey')
    
        # sns.scatterplot(data = df_using[df_using.index == 'pos1 C-G; pos2 A-T; pos4 G-C; pos5 T-A; '], # 0.460, 0.660
        #                 color  ='red', 
        #                 x = x, y = y, s = 15)
        plt.xlim(-0.05, 1)
        plt.ylim(-13, 4)
        plt.xticks([0, 0.25, 0.5, 0.75, 1], color = 'white', size =0) 
        plt.yticks([-12, -8, -4, 0, 4], color = 'white', size =0) 
        
        sns.despine()
        
        plt.axhline(cutoff_lce, color = 'k', lw = 1, ls = '--' )
        plt.axvline(cutoff_rc0, color = 'k', lw = 1, ls = '--' )
    
        plt.ylabel('',fontsize =0, color = 'white')
        plt.xlabel('',fontsize =0, color = 'white')
        
        plt.savefig(f'scatterplot lce0 vs rc0 all 6-bp {str(len(df_using))} variants.png', bbox_inches='tight', dpi =300)
        
    
    # calculate enrichment
    def calculate_enrichment():
        df_top = df_using[(df_using['rc0_dro'] >= cutoff_rc0) & (df_using['lce0_dro'] >= cutoff_lce) ].copy()
        
        df_enrichment = df_comb_6bp[df_comb_6bp['number_variant_dro'] > 0 ][['comb', 'number_variant_dro', 
                                                                             'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'group']].copy()
        df_enrichment = df_enrichment.set_index('comb')
        df_enrichment.columns = ['count_all', 'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'group', ]
        df_enrichment = df_enrichment[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'group', 'count_all']]
        
        df_enrichment['count_top'] = df_enrichment.index.map(lambda x: 
                                                             len(us.extract_variant_given_motif(df_top, x)))
    
        df_enrichment = df_enrichment[df_enrichment['count_top'] > 0]
        df_enrichment['group'].value_counts()
        
        df_enrichment['fraction_all'] = df_enrichment['count_all']/len(df_using)
        df_enrichment['fraction_top'] = df_enrichment['count_top']/len(df_top)
        df_enrichment['enrichment'] = df_enrichment.apply(lambda x: np.log2(x['fraction_top']/x['fraction_all']),
                                                          axis = 1 )
        df_enrichment.to_csv('./processed table/topmotif_enrichment.txt', sep = '\t', index = True)

    def draw_motif_enrichment_and_shifting_score():
        df_draw = pd.read_csv('./processed table/topmotif_enrichment.txt', 
                                     sep = '\t', index_col = 0).fillna('')
        
        df_draw = df_draw[df_draw['count_top'] >= 3]
        df_draw['shift_comb_count'] = df_shift_dro['shift_comb_count']
        df_draw['shift_score'] = df_shift_dro['shift_score']
        df_draw['rc0_dro'] = df_shift_dro['rc0_dro']
        
        cutoff_enrichment = 2
        cutoff_shift = 0.6
        
        # df_draw = df_draw[df_draw['shift_score'] >= 0.6]
        df_draw = df_draw[df_draw['shift_comb_count'] > 1]
        # df_draw = df_draw[df_draw['enrichment'] >= 1]
        
        df_draw['group'].value_counts()
        
        # df_draw = df_draw[df_draw['rc0_dro'] >= 0.5]
        
        fig, ax = plt.subplots(1, 1, figsize=(4,4))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
        ct.decorate_boxplot(ax)

        sns.scatterplot(data = df_draw.sort_values('group', ascending = False), hue = 'group',
                        x = 'enrichment', y = 'shift_score',  linewidth = 0.1, ec = 'k',
                        alpha = 1, s=15, palette = 'RdYlGn_r')


        plt.axvline(cutoff_enrichment, color = 'k', lw = 1, ls = '--' )
        plt.axhline(cutoff_shift, color = 'k', lw = 1, ls = '--' )
    
        plt.ylabel('',fontsize =0, color = 'white')
        plt.xlabel('',fontsize =0, color = 'white')
    
        plt.xticks([-6, -4, -2, 0, 2, 4], color = 'white', size =0) 
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2], color = 'white', size =0) 
        
        sns.despine()
    
        plt.legend('',frameon=False)

        plt.savefig('scatterplot shifting score vs motif enrichment.png', bbox_inches='tight', dpi =300)

    def swarmplot_dres_score():
        df_dres = us.import_dres()
        
        fig, ax = plt.subplots(1, 1, figsize=(4,4))
        sns.set_style('ticks', {'axes.edgecolor': 'black',
                                'grid.linestyle': 'dashed',
                                'axes.grid': True,})
        ct.decorate_boxplot(ax)

        sns.swarmplot(data = df_dres, x = 'group', y = 'dres_score',  s=5,
                      order = ['3 motif', '4 motif', '5 motif'], ec = 'k', linewidth = 0.2, 
                      palette=['#fffebe', '#b7e075', '#4bb05c'])
        plt.ylabel('',fontsize =0, color = 'white')
        plt.xlabel('',fontsize =0, color = 'white')
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)
        plt.yticks([0.4, 0.6, 0.8, 1, 1.2], color = 'white', size =0) 
        plt.xticks( color = 'white', size =0) 
        sns.despine()
    
        plt.legend('',frameon=False)

        plt.savefig('swarmplot dres score of 333 DRES.png', bbox_inches='tight', dpi =300)

#%% (DONE) FIGURE 4: SCATTER PLOT - SHIFT SCORE 
    
def compare_5bpdres_vs_6bpdres():
    df_using = df_all_6bp.copy()[['pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'rc0_dro', 'lce0_dro']].dropna()
    df_dres = us.import_dres()

    list_3bp_dres = df_dres[df_dres['group'] == '3 motif'].index.tolist()
    list_4bp_dres = df_dres[df_dres['group'] == '4 motif'].index.tolist()
    list_5bp_dres = df_dres[df_dres['group'] == '5 motif'].index.tolist()

    df_all_dres = pd.DataFrame()

    for motif in list_3bp_dres + list_4bp_dres +  list_5bp_dres:
        df =  us.extract_variant_given_motif(df_using, motif)
        
        # df['motif'] = motif
        df['position'] = ''.join([s for s in motif if s.isdigit()])
        df_all_dres = pd.concat([df, df_all_dres], axis = 0)

    df_all_dres= df_all_dres.reset_index().groupby(['index', 'pos1', 'pos2', 'pos3',
                                      'pos4', 'pos5', 'pos6', 'rc0_dro', 'lce0_dro'])['position'].apply(lambda x: x.sum()).reset_index()

    df_all_dres['position'] = df_all_dres['position'].map(lambda x: ''.join(sorted(list(set([i for i in x])))) )
    df_all_dres['len'] = df_all_dres['position'].map(lambda x: '6' if str(len(x)) == '6' else 'other' )
    
    df_all_dres['len'].value_counts()
    
    
    # draw 6-bp DRES
    value = 'rc0_dro'
    value = 'lce0_dro'
    
    fig, ax = plt.subplots(1, 1, figsize=(2,4))
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})
    boxprops = {'edgecolor': 'k', 'linewidth': 1.5}
    lineprops = {'color': 'k', 'linewidth': 1.5}
    boxplot_kwargs = {'boxprops': boxprops,  'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,}
    

    sns.boxplot(data = df_all_dres, x = 'len', y = value, order = ['other', '6'],
                fliersize=0.5, showfliers=True, **boxplot_kwargs,
                                palette=['#60d394', '#25a18e' ], saturation=1 )


    ct.adjust_box_widths(fig, 0.8)

    plt.xticks( fontsize =0, color = 'white')
    plt.xlabel('',fontsize =0, color = 'white')

    if value == 'rc0_dro':
        plt.ylim(-0.05,1)
        plt.yticks([0,0.2, 0.4, 0.6, 0.8, 1.0], fontsize =0, color = 'white') 

    else:
        plt.ylim(-10,4)
        plt.yticks([-8, -4, 0, 4], fontsize =0, color = 'white') 
        
    plt.ylabel('',fontsize =0, color = 'white')
    sns.despine()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ct.decorate_boxplot(ax)
    
    plt.savefig(f'boxplot compare 6-bp dres vs other dres {value}', bbox_inches= 'tight', dpi = 300)    

    # statistic
    value = 'rc0_dro' # 4.081432520027402e-86
    value = 'lce0_dro' # 7.702154328320636e-14
    
    ct.statistic_mannwhitneyu (df_all_dres[df_all_dres['len'] == 'other'][value],
                               df_all_dres[df_all_dres['len'] == '6'][value],
                               type_ = 'two-sided') 
    

def pieplot_fraction_bp_in_dres():
    df_using= df_shift_dro.copy()
    df_using = df_using[df_using['shift_comb_count'] != 1]
    df_using['group'].value_counts()

    number_bp = 5
    df_using = df_using[df_using['group'] == f'{number_bp} motif']
    

    fig, ax = plt.subplots(1, 1, figsize=(4,4))
    sns.set_style('ticks', {'axes.edgecolor': 'black',
                            'grid.linestyle': 'dashed',
                            'axes.grid': True,})

    
    df_top = us.import_dres()
    df_top = df_top[df_top['group'] == f'{number_bp} motif']
    df_top['group'].value_counts()

    df_enrichment = pd.DataFrame(index = ['A-T', 'T-A', 'C-G', 'G-C', 'G-T', 'T-G'])
    
    for i in range(1,7):
        df_enrichment[f'pos{i}_all'] = df_using[f'pos{i}'].value_counts()
        df_enrichment[f'pos{i}_top'] = df_top[f'pos{i}'].value_counts()
    df_enrichment['color'] = df_enrichment.index.map(lambda x: '#60d394' if x in ['C-G', 'G-C'] 
                                                    else ( '#3da4dc' if x in ['T-A', 'A-T'] else '#d3d3d3' ) )

    df_enrichment= df_enrichment.fillna(0)
    # pie plot
    p = 'pos1_top'
    for i in range(1,7):
        fig, ax = plt.subplots(1, 1, figsize=(4,4))
        p = f'pos{i}_top'
        plt.pie(df_enrichment[[p, 'color']].dropna().sort_values(p, ascending = False)[p].tolist() ,
                colors = df_enrichment[[p, 'color']].dropna().sort_values(p, ascending = False)['color'].tolist())
        plt.savefig(f"pieplot {p} {number_bp}-bp dres in {df_enrichment[p].sum()} motif.png",
                    bbox_inches="tight",dpi = 300)
