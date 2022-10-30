# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:52:50 2019

@author: Trung Duc Nguyen
"""

#%% IMPORT LIBRARY

import sys
import pandas as pd
import RNA
import forgi.utilities.stuff as fus
import forgi.graph.bulge_graph as cgb
import re

start = pd.Timestamp.now()

#%% 1) FOLD RNA USING RNAfold

df_input = pd.read_csv(sys.argv[1], sep = '\t', names = ['variant', 'seq'] ).set_index('variant')

df_output = df_input.copy()

def fold(seq):
    md = RNA.md()
    md.noLP = True
    fc = RNA.fold_compound(seq, md)
    return(fc.mfe()[0])


df_output['dot_struct'] = df_output['seq'].map(lambda x: fold(x))

print('Finish RNAfold...')


#%% 2) FORGI

def make_define_struct(dot_struct):
    x = cgb.BulgeGraph.from_dotbracket(dot_struct)
    return(x.to_element_string())

df_output['define_struct'] = df_output['dot_struct'].map(lambda x: make_define_struct(x))
df_output['forgi_struct'] = df_output['dot_struct'].map(lambda x: ' '.join([str(i) for i in fus.dotbracket_to_pairtable(x)]))

print('Finish forgi...')


#%% 3) MAKE NEW DEFINE STRUCTURE

'''
need to define:
    L: loop (h)
    M: match (s)
    B: single bulge // bulge in one strand (i)
    S: symetric bulge (i)
    A: asymetric bulge (i)
    U: multiple loop (m)
    F: 5p fragment
    T: 3p fragment

make define count from 1
     find the position of Mi{...}M
     run_start: position of first M, counted from 1
     run_end: position of last M, counted from 1
     from run_start, run_end: take look back the forgi list (forgi list containing which 
     position base pair with which position)
    
     calculate the distance from run_start to run_end --> the number of i at sense = distance sense - 1
     calculate the distance antisense --> the number of i at antisense = distance antisense - 1
     if distance antisense == distance sense --> same mismatch at both strand --> i at sense is S
     if distance_antisense != distance_sense and distance_antisense == 1 --> no mismatch at antisense --> i at sense is B
     if distance_antisense != distance_sense and distance_antisense != 1 --> different mismatches --> i at sense is A
'''

def make_list_i_new_define (forgi_struct, define_struct):
    list_i1 = []
    list_i2 = []
 
    list_start = [m.start() for m in re.finditer(r"Mi{1,100}", '0' + define_struct)] 
    list_end = [m.end() for m in re.finditer(r"i{1,100}M", '0' + define_struct)] 
    
    
    forgi_list = forgi_struct.split(' ')[:-1]
    forgi_list = [int(i) for i in forgi_list]

    for i in range(0,len(list_start)):
        run_start = list_start[i]
        run_end = list_end[i] - 1
        #print(str(run_start) + " to " + str(run_end))
        #print('base pair are:' + str(forgi_list[run_start]) + " to " + str(forgi_list[run_end]))
    
        distance_sense  = abs(run_start - run_end)
        distance_antisense = abs(forgi_list[run_start] - forgi_list[run_end])
        
        if distance_antisense == distance_sense:
            list_i1.append( (distance_sense -1) * 'S')
            list_i2.append( (distance_sense -1) * 'S')
        elif distance_antisense != distance_sense and distance_antisense == 1:
            list_i1.append( (distance_sense -1) * 'B')
            list_i2.append( (distance_sense -1) * 'B')    
        elif distance_antisense > distance_sense and distance_antisense != 1:
            list_i1.append( (distance_sense -1) * 'A')
            list_i2.append( (distance_sense -1) * 'S')
        elif distance_antisense < distance_sense and distance_antisense != 1:
            list_i1.append( (distance_sense -1) * 'A')
            list_i2.append( (distance_antisense -1) * 'S' + (distance_sense - distance_antisense) * 'A' )
            
    return(list_i1, list_i2)

'''
    define_structure replace i{1,100} to space and splitted to list
    add each element from list_new_define with each element from list_i
'''

def make_new_define_struct(forgi_struct, define_struct):
    replace_dict = {'f': 'F','t': 'T', 'h': 'L', 's': 'M', 'm': 'U'}
    for old, new in replace_dict.items():
        define_struct = define_struct.replace(old, new)
    
    list_i1, list_i2 = make_list_i_new_define(forgi_struct, define_struct)
    
    list_new_define_struct = re.sub(r'i{1,100}', ' ', define_struct).split(' ')
    new_define_struct1 = ''
    new_define_struct2 = ''
    
    for i in range(0, len(list_i1)):
        new_define_struct1 = new_define_struct1 + list_new_define_struct[i] + list_i1[i] 
        new_define_struct2 = new_define_struct2 + list_new_define_struct[i] + list_i2[i] 
    
    new_define_struct1 = new_define_struct1 + list_new_define_struct[-1]
    new_define_struct2 = new_define_struct2 + list_new_define_struct[-1]
    
    return(new_define_struct1, new_define_struct2)



# seq ='CAUUAUUACUUUUGGUACGCGCUGUGACACUUCAAACUCGUACCGUGAGUAAUAAUGCG'
# dot_struct = fold(seq)
# define_struct = make_define_struct(dot_struct)
# forgi_struct = ' '.join([str(i) for i in fus.dotbracket_to_pairtable(dot_struct)])
# new_define_struct1, new_define_struct2 = make_new_define_struct (forgi_struct, define_struct)


df_output['define_struct'] = df_output['dot_struct'].map(lambda x: make_define_struct(x))
df_output['forgi_struct'] = df_output['dot_struct'].map(lambda x: ' '.join([str(i) for i in fus.dotbracket_to_pairtable(x)]))

df_output[['new_define_struct1','new_define_struct2']] = df_output.apply(lambda x: make_new_define_struct(x['forgi_struct'], 
                                                                                                           x['define_struct']),
                                                                         result_type ='expand',
                                                                         axis = 1)

print('Finish making new structure...') 

#%% 4) MAKE NEW DEFINE STRUCTURE WOBBLE

def make_new_define_struct_with_wobble (seq, forgi_struct, new_define_struct2):

    seq = '0' + seq
    new_define_struct2 = '0' + new_define_struct2

    forgi_list = [int(i) for i in forgi_struct.split(' ')]

    new_define_struct_list = [i for i in new_define_struct2]

    new_define_struct_wobble_list = []

    # print(seq)
    for position in range(1, len(new_define_struct_list)):
        positional_struct = new_define_struct_list[position]
        if positional_struct != 'M':
            new_define_struct_wobble_list.append(positional_struct)
        elif seq[position] + seq[forgi_list[position]] in ['GT', 'GU'] :
            new_define_struct_wobble_list.append('W')
        elif seq[position] + seq[forgi_list[position]] in ['TG', 'UG'] :
            new_define_struct_wobble_list.append('w')
        else:
            new_define_struct_wobble_list.append('M')

    return( ''.join(new_define_struct_wobble_list))

df_output['new_define_struct_wobble'] = df_output.apply(lambda x: make_new_define_struct_with_wobble(x['seq'], 
                                                                                                   x['forgi_struct'], 
                                                                                                   x['new_define_struct2']), 
                                                      axis = 1)

print('Finish making new structure with wobbles...')  
 

#%% 5) CONCRETE STRUCTURE

def find_mismatch(struct):


    '''
    input: 
        struct = 'MMAAAAMMMMMMMSMMMMM'

    A, S, M, B, F, T, H

    output:
        position = 3, size = 4, type = AAAA 
        position = 14, size = 1, type = S

    --> return: 3AAA14S
    number show the starting position
    '''

    output = ''
    count = 0
    if struct.replace('M', '').replace('U', '').replace('T', '').replace('F', '')   == '':
        return(['none',0])
    else:
        for m in re.finditer(r"[A|S|B|L]{1,100}", ('0'+ struct).replace('T', '').replace('F', '')  ):
            start_pos = m.start()
            end_pos = m.end()
            types = ('0'+ struct).replace('T', '').replace('F', '')[start_pos:end_pos]
            output= output+str(start_pos).zfill(2) +'-'+types +' '
            count +=1
        return([output[:-1], count])


def reverse_mismatch(mismatch):
    '''
    for example: 04-S --> S-04
    '''
    reverse_mismatch =  mismatch.split('-')[1] + '-' + mismatch.split('-')[0]
    return(reverse_mismatch)

def find_position_bulge_another_strand (bulge, bulge_in, strand_5p, strand_3p):

    '''
    04-S S-04; 21-B-00; 26-B-00; 28-LLLL-33; 00-BBBBBBB-23
    04-S S-04; 21-B-(20); 26-B-(31); 28-LLLL-33; (20)-BBBBBBB-23

    strand_5p = 'FMMMSMMMMMMMMMMMMMMMMBMMMMBMLLLL'  
    
    strand_3p = 'TTTMMMSMMMMMMMMMMMMMMMMMMBBBBBBBMMMLLLL'
    bulge_in = '3p'
    bulge_in = '5p'
    bulge = '26-B-00'
    '''
    
    if bulge_in ==  '5p':
        count_M_till_bulge = strand_5p.replace('F', '')[:int(bulge.split('-')[0])].count('M')
        #print(strand_5p.replace('F', '')[:int(bulge.split('-')[0])])
        count_M_in_3p = 0
        count_nu_in_3p = 0
        for nu in strand_3p.replace('T', ''):
            if nu == 'M':
                count_M_in_3p += 1
            count_nu_in_3p +=1
            if count_M_in_3p == count_M_till_bulge:
                break
        return('{}-({})'.format(bulge,str(count_nu_in_3p).zfill(2)))
        
    if bulge_in ==  '3p':
        count_M_till_bulge = strand_3p.replace('T', '')[:int(bulge.split('-')[1])].count('M')
        #print(strand_3p.replace('T', '')[:int(bulge.split('-')[1])])
        count_M_in_5p = 0
        count_nu_in_5p = 0
        for nu in strand_5p.replace('F', ''):
            if nu == 'M':
                count_M_in_5p += 1
            count_nu_in_5p +=1
            if count_M_in_5p == count_M_till_bulge:
                break
        return('({})-{}'.format(str(count_nu_in_5p).zfill(2), bulge))

def join_5p_struct_and_3p_struct (struct_5p, struct_3p, strand_5p, strand_3p):
    '''
        strand_5p = 'FMMMSMMMSMMMMMMMMMMMMMBBBBBBBBBMMMLLLLLL'  
        strand_3p = 'TTTMMMSMMMSMMMMMMMMMMMMMMMMLLLLLL'
        
        struct_5p = '04-S 08-S 22-BBBBBBBBB 34-LLLLLL'
        struct_3p = '04-S 08-S 25-LLLLLL'
        struct_out = 04-S S-04; 08-S S-08; 22-BBBBBBBBB; 34-LLLLLL-25    
    
    '''
    
    list_struct_5p = struct_5p.split(' ')
    list_struct_3p = struct_3p.split(' ')
    list_out = []
    for mm5p in list_struct_5p:
        if 'S' in mm5p:
            for mm3p in list_struct_3p:
                if 'S' in mm3p:
                    list_out.append(mm5p + ' ' +  reverse_mismatch(mm3p))
                    list_struct_3p.pop(list_struct_3p.index(mm3p))
                    break
        elif 'A' in mm5p:
            for mm3p in list_struct_3p:
                if 'A' in mm3p:
                    list_out.append(mm5p + ' ' +  reverse_mismatch(mm3p))
                    list_struct_3p.pop(list_struct_3p.index(mm3p))
                    break
        elif 'L' in mm5p:
            for mm3p in list_struct_3p:
                if 'L' in mm3p:
                    list_out.append(mm5p + '-' + mm3p.split('-')[0])
                    list_struct_3p.pop(list_struct_3p.index(mm3p))
                    break
        elif 'B' in mm5p:
            list_out.append(find_position_bulge_another_strand (mm5p, '5p', strand_5p, strand_3p))
    for mm3p in list_struct_3p:
        if 'B' in mm3p:
            list_out.append(find_position_bulge_another_strand (reverse_mismatch(mm3p), '3p', strand_5p, strand_3p))
    return('; '.join(list_out))



def make_concrete_struct(new_define_struct2):
    '''
    new_define_structure must:
        have only one loop
    '''
#    new_define_structure =  'FFFFFFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMLLLMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMTTTTTTT'
    
    strand_5p = new_define_struct2.split('L')[0] + 'L'* new_define_struct2.count('L')
    strand_3p = new_define_struct2.split('L')[-1][::-1] + 'L'*new_define_struct2.count('L')
    
    struct_5p = find_mismatch(strand_5p)[0]
    struct_3p = find_mismatch(strand_3p)[0]
    
    struct_joined = join_5p_struct_and_3p_struct(struct_5p, struct_3p, strand_5p, strand_3p)
    struct_joined = str(strand_5p.count('F')) + 'F; ' + struct_joined + '; ' + str(strand_3p.count('T')) + 'T'
    return(struct_joined)

df_output['concrete_struct'] = df_output.apply(lambda x: 
                                               make_concrete_struct(x['new_define_struct2']) if x['new_define_struct2'].count('ML') == 1 else 'multiple loop', 
                                               axis = 1)

print('Finish making concrete structure...') 

#%% 6) SAVE RNA STRUCTURE

# print(df_output.head())
df_output.reset_index().to_csv( sys.argv[2], sep = '\t', header = True, index = False)

print(pd.Timestamp.now() - start)