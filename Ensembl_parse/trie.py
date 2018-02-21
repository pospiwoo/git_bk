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
 
class SuffixTrie:
    def __init__(self, prmtr):
        self.prmtr = prmtr
        self.root = TreeNode('s', 's')
        self.tree_node_len = 0
        
    def makeST(self, tree, seq, coor_list):
        if len(seq) != len(coor_list):
            print ('seq and coor_list length is different')
            return False
        for i in range(0,len(seq)-self.prmtr.max_idx):
            tree.submakeST(seq[i:i+self.prmtr.max_idx], coor_list[i:i+self.prmtr.max_idx])
                                
    def submakeST(self, seq, coor_list):
        cur = self.root
        for i in range(0, len(seq)):
            if seq[i] not in cur.out:
                for j in range(i, len(seq)-1):
                    cur.out[seq[j]] = TreeNode(seq[j], coor_list[0])
                    cur = cur.out[seq[j]]
                j = len(seq)-1
                cur.out[seq[j]] = TreeNode(seq[j], coor_list[0])
                break
            else:
                cur.coor[coor_list[0]] = True
                cur = cur.out[seq[i]]
                continue

    def searchPath(self, seq):
        cur = self.root
        for i in range(0, len(seq)-1):
            if seq[i] in cur.out:
                cur = cur.out[seq[i]]
                continue
            else:
                return False
        if seq[-1] in cur.out:
            return cur.out[seq[-1]]
                
    def retrievePath(self, cur_node, cur_seq):
        if cur_node.out == {}:
            self.printTreeFile.write( cur_seq + ' -> ' + cur_node.lab + '\n' )
            return True
        for idx in cur_node.out:
            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab)
#            self.retrievePath(cur_node.out[idx], cur_seq + cur_node.lab + ' -> ')
                    
    def printTree(self):
        cur = self.root
        with open('F:\\Ensembl\\printTree.txt','w') as self.printTreeFile:
            self.retrievePath(cur, '')
    
        