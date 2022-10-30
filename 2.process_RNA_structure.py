# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:15:45 2022

@author: Trung Duc Nguyen
"""
#%% (DONE) IMPORT LIBRARY

import pandas as pd
import common_tools as ct
import using_now as us

#%% (DONE) PROCESS RNA STRUCTURE

def find_main_loop(split_loop_structure):
    list_stem_loop = split_loop_structure.split(' ')
    main_loop =''
    for i in list_stem_loop:
        if 'L' in i:
            count_m = i.count('M')/2
            if count_m > main_loop.count('M')/2:
                main_loop = i
    return(main_loop)
    
def concrete_main_loop(main_loop):
    concreted_main_loop = '{}-[{}]-{}'.format(
            main_loop.split('-')[0],
            ct.concrete_structure(main_loop.split('-')[1]),
            main_loop.split('-')[2])
    return(concreted_main_loop)

def reconcrete_struct():

    df_struct = pd.read_csv('./processed_table/576_structure.txt', sep = '\t')
    
    df_struct['Split_multiple_loop_structure'] = df_struct.apply(lambda x: 
        ct.split_multiple_loop(x['Sequence'], x['Dot_structure'], x['New_define_structure']) , axis = 1)
    
    df_struct['Main loop'] = df_struct['Split_multiple_loop_structure'].map(lambda x: find_main_loop(x))
        
    df_struct['Main loop concreted'] = df_struct['Main loop'].map(lambda x: concrete_main_loop(x))

    df_struct = df_struct.set_index('Variant')
    
    df_struct[['Sequence', 'New_define_structure', 'New_define_structure_wobble', 
               'Split_multiple_loop_structure', 'Main loop', 'Main loop concreted']].to_csv('./processed_table/576_structure_concrete.txt', sep = '\t')


def process_struct_with_classification():
    df_struct = us.import_struct()[['Main loop concreted', 'Name_struct']]
    
    df_sg = us.import_sg()  # second type of subgroup
    df_struct = pd.concat([df_struct, df_sg], axis = 1) 
    
    df_struct.columns = ['main_loop', 'name_struct', 'ref_seq', 'sg', 'motif']

    def classify_struct(main_loop):
        ''' 
        - if main loop of structure contains '19-' and '-90' --> good struct
        - among good struct: 
            if no 'B' no 'A' --> symm struct
            else --> asymm
        '''
        if ('19-' in main_loop) and ('-90' in main_loop):
            if ('A' not in main_loop) and ('B' not in main_loop):
                return('symm')
            else:    
                return('asymm')
        else:
            return('unknown')

    df_struct['type_struct'] = df_struct['main_loop'].map(lambda x: classify_struct(x))
    
    '''
    short_name: remove other common structure
    short_name_from_rand: count the mismatch from the randomized position
    '''
    
    df_struct['short_name'] = df_struct.apply(lambda x: x['main_loop'].replace('17-S S-17; 30-S S-30; 34-LLLLLL-34', '').replace('19-[', '').replace(']-90', '').replace(' ;', '') if 
                                              (x['type_struct'] != 'unknown') else '', axis = 1 )
    
    df_struct['short_name_from_rand'] = df_struct['short_name'].map(lambda x: 
                                            x.replace('3', '1').replace('4', '2').replace('5', '3').replace('6', '4').replace('7', '5').replace('8', '6').replace('9', '7'))

    df_struct['pos1'] = df_struct['motif'].map(lambda x: x[0]+ '-' +x[12])
    df_struct['pos2'] = df_struct['motif'].map(lambda x: x[1]+ '-' +x[11])
    df_struct['pos3'] = df_struct['motif'].map(lambda x: x[2]+ '-' +x[10])
    df_struct['pos4'] = df_struct['motif'].map(lambda x: x[3]+ '-' +x[9])
    df_struct['pos5'] = df_struct['motif'].map(lambda x: x[4]+ '-' +x[8])
    df_struct['pos6'] = df_struct['motif'].map(lambda x: x[5]+ '-' +x[7])
    
    df_struct.to_csv('processed_table/576_struct_with_classification.txt', sep ='\t')
    return()