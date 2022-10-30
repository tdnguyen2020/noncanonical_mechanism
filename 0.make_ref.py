# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:15:40 2021

@author: Trung Duc Nguyen
"""

#%%                

import itertools
import pandas as pd

#%%

def generate_dict_backbone():
    list_all_6N = [''.join((i)) for i in list(itertools.product(['R', 'Y'], repeat=6))]
    
    dict_backbone = {}
    
    for i in list_all_6N:
        
        backbone = '{}-{}'.format(i, i.replace('R', 'y').replace('Y', 'r').upper()) 
        list_individual = []
        for nu in backbone:
            if nu == 'R':
                list_individual.append(['A', 'G'])
            elif nu == 'Y':
                list_individual.append(['C', 'T'])
            else:
                list_individual.append(' ')
        
        list_variant_subgroup = []
        for nu in itertools.product(*list_individual):
            list_variant_subgroup.append('CTGCCATTTTACAATCCAAC{}TCTAATTTCTCCACGTCTTTGGTAATAAGGTTTGGCAAAGATGTGGAAAAATTGGA{}GTCTTTCAACCCAACCGGAA'.format(''.join(nu)[:6], ''.join(nu)[-6:][::-1]))
    #    print(list_individual)
        dict_backbone['{}'.format(i)] = list_variant_subgroup
    return(dict_backbone)
    
def make_ref():
    dict_backbone = generate_dict_backbone()
    
    list_all_variant = []
    for value in dict_backbone.values():
        list_all_variant = list_all_variant + value
    
    file_fasta_output = open('./processed_table/576_ref.fa', 'w')
    
    for variant in list_all_variant:
        file_fasta_output.write('>Variant_' + str(list_all_variant.index(variant)+1).zfill(6) + '\n')
        file_fasta_output.write(variant + '\n')
    
    file_fasta_output.close()
    file_txt_output = open('./processed_table/576_ref.txt', 'w')
    
    for variant in list_all_variant:
        file_txt_output.write('Variant_' + str(list_all_variant.index(variant)+1).zfill(6) + 
                              '\t' + variant + '\n')
    file_txt_output.close()


def make_allSG():
    dict_backbone = generate_dict_backbone()
    
    df_refseq = pd.read_csv('./processed_table/576_ref.txt', sep ='\t', names = ['Variant', 'Ref_seq']).set_index('Ref_seq')
    
    for key, value in dict_backbone.items():
        df_refseq[key] = df_refseq.index.map(lambda x: 
            x.replace('CTGCCATTTTACAATCCAAC', '').replace('GTCTTTCAACCCAACCGGAA', '').replace('TCTAATTTCTCCACGTCTTTGGTAATAAGGTTTGGCAAAGATGTGGAAAAATTGGA', '-') if x in dict_backbone[key] else 'no')
    df_refseq = df_refseq.reset_index().set_index('Variant')
    df_refseq.to_csv( './processed_table/576_ref_allSG.txt', sep ='\t')  


def make_longtable():
    df_refseq = pd.read_csv('./processed_table/576_ref.txt', sep ='\t', names = ['Variant', 'Ref_seq']).set_index('Ref_seq')
    df_refseq['randomized']  = df_refseq.index.map(lambda x: 
            x.replace('CTGCCATTTTACAATCCAAC', '').replace('GTCTTTCAACCCAACCGGAA', '').replace('TCTAATTTCTCCACGTCTTTGGTAATAAGGTTTGGCAAAGATGTGGAAAAATTGGA', '-'))
        
    df_refseq.to_csv( './processed_table/576_ref_allSG_longtable.txt', sep ='\t', index = False) 


def make_allSG_2nd():
    df_sg = pd.read_csv('./processed_table/576_ref_allSG.txt', sep ='\t').set_index('var')
    df_sg = df_sg.reset_index().melt(id_vars=['var', 'ref_seq'])
    df_sg = df_sg[df_sg['value'] != 'no']
    df_sg.columns = ['var', 'ref_seq', 'sg', 'motif']
    df_sg = df_sg.set_index('var')
    df_sg.to_csv('./processed_table/576_ref_allSG_2nd.txt',  sep ='\t')   

