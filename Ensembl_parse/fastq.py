# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 20:42:47 2017

@author: Sunghee Woo
"""
#import os, sys, pickle
import parameter, gene
import random
from collections import defaultdict


class ProcessFASTQ(object):
    def __init__(self):
        self.found = 0
        self.notfound = 0
        self.multiple_results_cnt = defaultdict(int)
        self.prmtr = parameter.Parameter()
#        self.graph = gene.ConstructGraph()
        self.graph = gene.unpickleGraph()
        self.parseFASTQ()
        
    def parseFASTQ(self):
        print ('processing FASTQ')
        with open(self.prmtr.fastq,'r') as in_fastq:
            for line in in_fastq:
                if (self.found + self.notfound) % 100000 == 0:
                    print ('processed:', self.found + self.notfound, 'found:',\
                           self.found, 'notfound:', self.notfound)
                    print (self.multiple_results_cnt)
                if line.startswith('#'):
                    continue
                if line.startswith('@'):
                    foundflag = False
                    seq = next(in_fastq).strip()
                    if self.prmtr.gen_rand_error:
                        seq = self.randSeqError(seq)
                    next(in_fastq) # info = 
                    next(in_fastq) # score = 
#                    found_flag1, coor1 =  self.graph.tree.searchPath(seq)
                    found_flag1, coor1 =  self.graph.tree.searchPath(seq[:50])
#                    found_flag2, coor2 =  self.graph.tree.searchPath(seq[50:])
##                    found_flag1, coor1 =  self.graph.tree.searchPath(seq[:31])
##                    found_flag2, coor2 =  self.graph.tree.searchPath(seq[31:62])
##                    found_flag3, coor3 =  self.graph.tree.searchPath(seq[62:])
##                    if found_flag1 and (not found_flag2):
##                        print (found_flag1, coor1, found_flag2, coor2)
##                    elif (not found_flag1) and found_flag2:
##                        print (found_flag1, coor1, found_flag2, coor2)
##                    elif found_flag1 and found_flag2:
##                        print (found_flag1, coor1, found_flag2, coor2)
##                    else:
##                        print (found_flag1, coor1, found_flag2, coor2)
                    if found_flag1:
#                        self.iterateGenome(seq[:50], coor1)
                        foundflag = self.iterateGenome(seq, coor1)
#                        self.iterateGenome(seq[:110], coor1)
#                    if found_flag2:
##                        self.iterateGenome(seq[:50], coor1)
#                        foundflag = self.iterateGenomeReverse(seq[:51], coor2)
##                        self.iterateGenome(seq[:110], coor1)
                    if foundflag:
                        self.found += 1
                    else:
                        self.notfound += 1
#                        print ('X', coor2, ':', seq[:50], seq[50:])
#                    if len(seq) < 60:
#                        continue
#                    if seq == '' or info == '' or score == '':
#                        continue
#                    found_flag3, coor3 =  self.graph.tree.searchPath(seq[:33])
#                    if found_flag3:
##                        self.iterateGenome(seq[:50], coor1)
##                        self.iterateGenome(seq, coor1)
#                        self.iterateGenome(seq[:50], coor3)
                    continue
                
    def iterateGenome(self, seq, coor):
        tmp = -1
        for i in coor:
            flag, returnlist = self.retrieveGenome(seq, i, 0)
            tmp += 1
            if flag:
#                self.found += 1
#                print (tmp, i, 'oooooo', returnlist)
#                print (tmp, i, 'oooooo')
                self.multiple_results_cnt[len(returnlist)] += 1
                return True
#            else:
#                self.notfound += 1
##                print (tmp, ':', i, 'xxxxxx', seq, returnlist)
        return False
                
    def retrieveGenome0(self, coordinate, err):
        if coordinate not in self.graph.GGraph.G:
            return False
#        cur_node = self.graph.GGraph.G[coordinate]
#        if cur_node.to_edges == {} or l_bases <= 1:
#            self.printTreeFile.write( cur_seq + cur_node.lab + '\n' )
#            return True
#        for idx in cur_node.out:
#            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
#            self.retrieveGenome()
            
    def retrieveGenome(self, seq, coordinate, err):
        return_list = []
        if coordinate not in self.graph.GGraph.G:
            return False, '1111111'
        cur_node = self.graph.GGraph.G[coordinate]
        if len(seq) == 1 and cur_node.base == seq[0]:
            return True, coordinate
        if len(seq) == 1 and cur_node.base != seq[0]:
            return False, []
        if cur_node.base != seq[0]:
#            err += 1
            return False, []
#        if err > self.prmtr.max_mismatch:
#            return False, '2222222'
        if cur_node.base == '':
#            err += 1
#            print ('ttttttttttt', coordinate, cur_node.base, seq)
            return False, '6546546464564'
        if cur_node.to_edges == {} and len(seq) > 1:
            return False, '33333333'
        for idx in cur_node.to_edges:
            flag, sub_coor = self.retrieveGenome(seq[1:], idx, err)
            if flag == False:
#                print (sub_coor)
                continue
            if sub_coor == []:
                continue
            if isinstance(sub_coor, int) or isinstance(sub_coor, str):
                tuple_list = [coordinate, sub_coor]
                return_list.append(tuple_list)
                continue
            else:
                if len(sub_coor) > 1:
                    for i in range(0,len(sub_coor)):
                        tuple_list = [coordinate]
                        for j in range(0,len(sub_coor[i])):
                            tuple_list.append(sub_coor[i][j])
                        return_list.append(tuple_list)
                else:
                    tuple_list = [coordinate]
                    for j in range(0,len(sub_coor[0])):
                        tuple_list.append(sub_coor[0][j])
                    return_list.append(tuple_list)
        if return_list == []:
            return False, '55555555'
        return True, return_list
            
    def retrieveGenome1(self, seq, coordinate, err):
        return_list = []
        if coordinate not in self.graph.GGraph.G:
            return False, '1111111'
        cur_node = self.graph.GGraph.G[coordinate]
        if len(seq) == 1 and cur_node.base == seq[0]:
            return True, coordinate
        if len(seq) == 1 and (err + 1 <= self.prmtr.max_mismatch):
            return True, coordinate
        if cur_node.base != seq[0]:
            err += 1
#            return False, []
        if err > self.prmtr.max_mismatch:
            return False, '2222222'
        if cur_node.to_edges == {} and len(seq) > 1:
            return False, '33333333'
        for idx in cur_node.to_edges:
            flag, sub_coor = self.retrieveGenome(seq[1:], idx, err)
            if flag == False:
#                print (sub_coor)
                continue
            if sub_coor == []:
                continue
            if isinstance(sub_coor, int) or isinstance(sub_coor, str):
                tuple_list = [coordinate, sub_coor]
                return_list.append(tuple_list)
                continue
            else:
                if len(sub_coor) > 1:
                    for i in range(0,len(sub_coor)):
                        tuple_list = [coordinate]
                        for j in range(0,len(sub_coor[i])):
                            tuple_list.append(sub_coor[i][j])
                        return_list.append(tuple_list)
                else:
                    tuple_list = [coordinate]
                    for j in range(0,len(sub_coor[0])):
                        tuple_list.append(sub_coor[0][j])
                    return_list.append(tuple_list)
        if return_list == []:
            return False, '55555555'
        return True, return_list
                
    def iterateGenomeReverse(self, seq, coor):
        tmp = -1
        for i in coor:
            flag, returnlist = self.retrieveGenomeReverse(seq, i, 0)
            tmp += 1
            if flag:
#                print (tmp, i, 'oooooo', returnlist)
#                print (tmp, i, 'oooooo')
                self.multiple_results_cnt[len(returnlist)] += 1
                return True
#            else:
#                print (tmp, ':', i, 'xxxxxx', seq, returnlist)
        return False
            
    def retrieveGenomeReverse(self, seq, coordinate, err):
        return_list = []
        if coordinate not in self.graph.GGraph.G:
            return False, '1111111'
        cur_node = self.graph.GGraph.G[coordinate]
#        print ('********', coordinate, cur_node.base, cur_node.from_edges)
        if len(seq) == 1 and cur_node.base == seq[-1]:
            return True, coordinate
        if len(seq) == 1 and cur_node.base != seq[-1]:
            return False, []
        if cur_node.base == '':
#            err += 1
#            print ('ttttttttttt', coordinate, cur_node.base, seq)
            return False, '6546546464564'
        if cur_node.base == 's':
#            err += 1
#            print ('ttttttttttt', coordinate, cur_node.base, seq)
            return False, '1212121212'
        if cur_node.base != seq[-1]:
#            err += 1
            return False, '9999999999999'
#        if err > self.prmtr.max_mismatch:
#            return False, '2222222'
        if cur_node.from_edges == {} and len(seq) > 1:
            return False, '33333333'
        for idx in cur_node.from_edges:
            flag, sub_coor = self.retrieveGenomeReverse(seq[:-1], idx, err)
            if flag == False:
#                print(1231365465431456456456456)
#                print (1231365465431456456456456, sub_coor)
                continue
            if sub_coor == []:
                continue
            if isinstance(sub_coor, int) or isinstance(sub_coor, str):
#                print(666666666666)
                tuple_list = [coordinate, sub_coor]
                return_list.append(tuple_list)
                continue
            else:
#                print(8888888888888)
                if len(sub_coor) > 1:
                    for i in range(0,len(sub_coor)):
                        tuple_list = [coordinate]
                        for j in range(0,len(sub_coor[i])):
                            tuple_list.append(sub_coor[i][j])
                        return_list.append(tuple_list)
                else:
                    tuple_list = [coordinate]
                    for j in range(0,len(sub_coor[0])):
                        tuple_list.append(sub_coor[0][j])
                    return_list.append(tuple_list)
        if return_list == []:
            return False, '55555555'
        return True, return_list
    
    def randSeqError(self, seq):
        new_seq = ''
        change_idx_list = random.sample(range(0, len(seq)), 2)
        for i in range(0,len(seq)):
            if i in change_idx_list:
                if seq[i] == 'A':
                    new_seq += 'T'
                elif seq[i] == 'T':
                    new_seq += 'C'
                elif seq[i] == 'C':
                    new_seq += 'G'
                elif seq[i] == 'G':
                    new_seq += 'A'
            else:
                new_seq += seq[i]
#        print (seq, new_seq)
        return new_seq
    
if __name__=='__main__':
    parse_obj = ProcessFASTQ()
    
    
    
    
    