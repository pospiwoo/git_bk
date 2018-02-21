# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 20:46:35 2017

@author: Sunghee Woo
"""
class TreeNode:
      def __init__(self, lab_str, idx):
           self.lab = lab_str
           self.coor = {}
           if idx != None:
               self.coor[idx] = True
           self.out = {}

class SuffixTree:
    def __init__(self, prmtr):
        self.prmtr = prmtr
        self.root = ''
        self.tree_node_len = 0
        
    def makeST(self, tree, seq, coor_list):
        if len(seq) != len(coor_list):
            print ('seq and coor_list length is different')
            return False
        if self.root == '':
            self.root = TreeNode('s', 's')
        for i in range(0,len(seq)-self.prmtr.max_idx):
#            if seq[i:i+self.prmtr.max_idx].startswith('CAAGAAGTGGAAGAAGACCAAGACCCATCATGCCCCAGGCTCAGCAGGGAGCTGCTGGAT'):
#                print ('', coor_list[i], seq[i:i+self.prmtr.max_idx])
            tree.submakeST(seq[i:i+self.prmtr.max_idx], coor_list[i:i+self.prmtr.max_idx])
        for i in range(len(seq)-self.prmtr.max_idx,len(seq)-self.prmtr.min_idx):
            tree.submakeST(seq[i:], coor_list[i:])
        
    def submakeST(self, seq, coor):
        cur = self.root
        i = 0
        while i < len(seq):
            if seq[i] not in cur.out:
                cur.out[seq[i]] = TreeNode(seq[i:], coor[0])
                return
            elif seq[i] in cur.out and len(cur.out[seq[i]].lab) == 1:
#                cur.coor[coor[0]] = True
                cur = cur.out[seq[i]]
                i += 1
                continue
            else:
                cur_child = cur.out[seq[i]]
                exists_flag = True
                j = 0
                while j < len(cur_child.lab) and j < len(seq)-i:
                    if seq[i+j] != cur_child.lab[j]:
                        exists_flag = False
                        pivot = j
                        break
                    j += 1
                if exists_flag:
                    cur_child.coor[coor[0]] = True
                    if len(cur_child.lab) > j and len(seq)-i > j:
                        cur = cur_child
                        i += j
                        continue
                    elif len(seq)-i == j:
                        return
                    elif len(seq)-i > j:
                        cur = cur_child
                        i += j
                        continue
                    else:
                        cur.coor[coor[0]] = True
                        print ('should not reach this point', cur_child.lab, seq, i, j)
                        return
                else:
                    if pivot == 0:
                        print ('wrong', seq, cur_child.lab, i, pivot)
                        return
                    end_node = TreeNode(seq[i+pivot:], coor[0]) # create end node (new sequence)
                    start_node = TreeNode(cur_child.lab[:pivot], coor[0]) # create start node (original sequence)
                    for idx in cur_child.coor:
                        if idx != 's':
                            start_node.coor[idx] = True
                    start_node.out[seq[i+pivot]] = end_node
                    start_node.out[cur_child.lab[pivot]] = cur_child
                    cur_child.lab = cur_child.lab[pivot:] # original childs label is curtailed
#                    cur_child.coor[idx] = True
                    cur.out[seq[i]] = start_node
                    break
    
    def searchPath(self, seq):
        cur = self.root
        i = 0
        while i < len(seq):
            if seq[i] not in cur.out:
                return False, []
            elif seq[i] in cur.out and len(cur.out[seq[i]].lab) == 1:
                cur = cur.out[seq[i]]
                i += 1
                continue
            else:
                cur_child = cur.out[seq[i]]
                j = 0
                while j < len(cur_child.lab) and j < len(seq)-i:
                    if seq[i+j] != cur_child.lab[j]:
                        return False, []
                    j += 1
                if len(cur_child.lab) > j and len(seq)-i > j:
                    cur = cur_child
                    i += j
                    continue
                elif len(seq)-i == j:
                    return True, cur_child.coor
                elif len(seq)-i > j:
                    cur = cur_child
                    i += j
                    continue
                else:
                    print ('should not reach this point', cur_child.lab, seq, i, j)
                    return False, []
        # check if last return is valid
        return True, cur_child.coor
    

    def retrievePath(self, cur_node, cur_seq):
        if cur_node.out == {}:
#            self.printTreeFile.write( cur_seq + cur_node.lab + '\n' )
            self.printTreeFile.write( cur_seq + cur_node.lab + ':' + ','.join(map(str, cur_node.coor)) + '\n' )
            return True
        for idx in cur_node.out:
#            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
            if cur_node.lab == 's':
                self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
            else:
                self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ':' + ','.join(map(str, cur_node.coor)) + ' -> ')

    def printTree(self):
        cur = self.root
        with open('F:\\Ensembl\\printTree.txt','w') as self.printTreeFile:
            self.retrievePath(cur, '')
            

#        print ('here', 0, seq[0:0+self.prmtr.max_idx])
#        tree.submakeST(seq[0:0+self.prmtr.max_idx], coor_list[0:0+self.prmtr.max_idx])
#        print ('here', 1, seq[1:1+self.prmtr.max_idx])
#        tree.submakeST(seq[1:1+self.prmtr.max_idx], coor_list[1:1+self.prmtr.max_idx])

#        seq += '$'
#        coor.append(-1)

#    def retrievePath(self, cur_node, cur_seq):
##        print (cur_node)
#        if cur_node.out == {}:
##            self.printTreeFile.write( cur_seq + ' s -> ' + cur_node.lab + '\n' )
##            if cur_seq == '' or cur_node.lab == '':
##                print (cur_node.lab)
##            self.printTreeFile.write( cur_seq + ' -> ' + cur_node.lab + '\n' )
##            print( ' -> ' + cur_node.lab + '\n' )
#            self.printTreeFile.write( cur_seq + cur_node.lab + '\n' )
#            return True
#        for idx in cur_node.out:
##            if idx == '' or cur_node.lab == {} : #or cur_node.lab == 's':
##                print (idx, cur_node.lab, cur_seq, cur_node.out)
##            print (idx, cur_node.lab)
#            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
##            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -' + idx + '> ')
##            if cur_node.lab == 's':
##                self.retrievePath(cur_node.out[idx], cur_seq)
##            else:
##                self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -' + idx + '> ')
##            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab)
#
#
#    def retrievePath1(self, cur_node, cur_seq):
#        if cur_node.out == {}:
#            self.printTreeFile.write( cur_seq + ' -> ' + cur_node.lab + '\n' )
#            return True
#        for idx in cur_node.out:
#            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab)
##            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
