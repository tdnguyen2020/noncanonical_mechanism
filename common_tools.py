# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 13:58:08 2020

@author: Trung Duc Nguyen


"""

#%%

import re
import forgi.graph.bulge_graph as fgb
import matplotlib.pyplot as plt 
import scipy.stats as ss
from statistics import mean, stdev, sqrt

#%%

def calculate_cohend(a, b):
    # test conditions
    # c0 = [2, 4, 7, 3, 7, 35, 8, 9]
    # c1 = [i * 2 for i in c0]
    
    cohens_d = (mean(a) - mean(b)) / (sqrt((stdev(a) ** 2 + stdev(b) ** 2) / 2))
    
    return(cohens_d)

def statistic_mannwhitneyu (array1,  array2, type_ = 'two-sided'):
    stat, p = ss.mannwhitneyu(array1, array2, alternative= type_)
    print(len(array1))
    print(len(array2))
    return (str(p))

def decorate_boxplot(ax):
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    plt.setp(ax.yaxis.get_ticklines(), 'markersize', 5)
    plt.setp(ax.yaxis.get_ticklines(), 'markeredgewidth', 2)
    plt.setp(ax.xaxis.get_ticklines(), 'markersize', 5)
    plt.setp(ax.xaxis.get_ticklines(), 'markeredgewidth', 2)
    plt.setp(ax.artists, edgecolor = 'k',)
    plt.setp(ax.lines, color='k')
    return()


def joined_2_lists_overlap(list_start, list_add):  
    
    '''
    given 2 lists.
    I want to 2 lists together if they share the same defined information
    '''
    
    list_overlap = [i for i in list_add if i in list_start ]        
    if list_overlap == []:
        return(list_start)
    else:
        return(list(set(list_start+ list_add)))


def split_multiple_loop (sequence, dot_structure, sequence_or_structure ):
    '''
    sequence_or_structure parameter:
        I can put any types: sequence, new_define_structure,...
    '''
            

    bg = fgb.BulgeGraph.from_dotbracket(dot_structure, sequence)


    list_connect_each_part = []
    dict_define_each_part  = {}

    # generate infromation: define each part of structure and connect them 
    
    defined_information =  bg.to_bg_string().split('\n')

    for line in defined_information:
        if 'define' in line:
            line =  line.replace('define ', '')
            line = line.split(' ')
            dict_define_each_part[line[0]] = line[1:]
        elif 'connect' in line:
            line = line.replace('connect ', '')
            line = line.split(' ')
            list_connect_each_part.append(line)

    # sort connect
    list_index = [int(i[0].replace('s', '')) for i in list_connect_each_part ]
    list_connect_each_part = [x for _,x in sorted(zip(list_index,list_connect_each_part))]

    # divide 
    list_position_stem_loop = []

    for i in list_connect_each_part:
        
        if 'h' in ''.join(i):
            stem_loop = i.copy()
        
            for j in list_connect_each_part[::-1]:
                stem_loop = joined_2_lists_overlap(stem_loop, j)
                if 'm' in ''.join(stem_loop):
                    break            
            # find smallest sx
#            print(stem_loop)
            min_sx = 's' + str(min(int(sx.replace('s', '')) for sx in stem_loop if 's' in sx)) 
#            print(min_sx)
            
            list_position_stem_loop.append( int( dict_define_each_part[min_sx][0])) 
            list_position_stem_loop.append( int( dict_define_each_part[min_sx][3]))

    list_position_stem_loop = sorted(list_position_stem_loop)
    list_position_stem_loop.append(len(sequence_or_structure)+1)
    
    
    output  = []
    
    if list_position_stem_loop[0] != 1:
        output.append('1-{}-{}'.format(sequence_or_structure[:list_position_stem_loop[0]-1], list_position_stem_loop[0]-1))
#    
    
    for i in range(len(list_position_stem_loop)-1):
        start = list_position_stem_loop[i]
        end = list_position_stem_loop[i+1]
    
        if i % 2 == 0:    
            output.append('{}-{}-{}'.format(str(start), sequence_or_structure[start-1:end], str(end)))
        else:
            output.append('{}-{}-{}'.format(str(start+1), sequence_or_structure[start:end-1], str(end-1)))
    
    return(' '.join([i for i in output if '--' not in i]))

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



def concrete_structure (new_define_structure):
    '''
    new_define_structure must:
        have only one loop
    '''
#    new_define_structure =  'FFFFFFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMLLLMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMTTTTTTT'
    
    strand_5p = new_define_structure.split('L')[0] + 'L'*new_define_structure.count('L')
    strand_3p = new_define_structure.split('L')[-1][::-1] + 'L'*new_define_structure.count('L')
    
    struct_5p = find_mismatch(strand_5p)[0]
    struct_3p = find_mismatch(strand_3p)[0]
    
    struct_joined = join_5p_struct_and_3p_struct(struct_5p, struct_3p, strand_5p, strand_3p)
    
    return(struct_joined)


#%%
if __name__ == '__main__':
    concrete_structure()

#%%
    
#import pandas as pd
#
#test = pd.read_csv('test_structure.txt.txt', sep = '\t')[['Variant', 'Sequence', 'Structure']].set_index('Variant')
#
#test['concrete structure'] = test['Structure'].map(lambda x: concrete_structure(x))
















