# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 20:45:30 2017

@author: Sunghee Woo
"""
from collections import defaultdict

class GraphNode:
      def __init__(self, base_str, coverage, node_type):
           #self.idx = idx
           self.base = base_str
           self.cov = coverage
           self.from_edges = defaultdict(int)
           self.to_edges = defaultdict(int)
           self.visited = 0
           self.type = node_type
           
class GGraph:
    def __init__(self):
        self.G = {}
        self.SNP = []
        self.SUB = []
        self.INS = []
        self.DEL = []
        self.graph_node_len = 0
    
    def addNode(self, idx, base_str, coverage, node_type):
        try:
            idx = int(idx)
        except ValueError:
            #print "Node idx should be INT", idx
            return False
        if base_str == 'N':
            print ('this base is masked', idx, base_str)
        if idx in self.G:
            #print "Node already exist, adding coverage to existing node", idx, base_str, coverage
            self.G[idx].cov += coverage
        else:
            self.G[idx] = GraphNode(base_str, coverage, node_type)
            self.graph_node_len += 1
        return True
    
    def addMutation(self, base_str, coverage, node_type, prev_ind, next_ind, ref_coor, cov):
        node_inst = GraphNode(base_str, coverage, node_type)
        if node_type == 'SNP':
            self.SNP.append(node_inst)
        elif node_type == 'SUB':
            self.SUB.append(node_inst)
        elif node_type == 'INS':
            self.INS.append(node_inst)
        elif node_type == 'DEL':
            self.DEL.append(node_inst)
        self.graph_node_len += 1
        self.addMutationEdge(node_inst, prev_ind, ref_coor, cov)
        self.addMutationEdge(node_inst, ref_coor, next_ind, cov)
        return True
    
    def addEdge(self, fro_idx, to_idx, cov):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ("Node indexes should be INT", fro_idx, to_idx)
            return False
        if fro_idx in self.G and to_idx in self.G:
            self.G[fro_idx].to_edges[to_idx] += cov # Assign edge from previous node
            self.G[to_idx].from_edges[fro_idx] += cov # Assign edge to next node
        else:
            print ("Nodes not exist", fro_idx, to_idx)
            return False
        return True
    
    def addMutationEdge(self, node_inst, fro_idx, to_idx, cov):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ("Node indexes should be INT", fro_idx, to_idx)
            return False
        node_inst.to_edges[to_idx] += cov # Assign edge from previous node
        node_inst.from_edges[fro_idx] += cov # Assign edge to next node
        return True
    
    def removeGraphNode(self, node_inst):
        node_inst.type += '_removed'
        for idx, cov in node_inst.from_edges.iteritems():
            from_node = self.G[idx]
            from_node.to_edges.pop(node_inst.idx, None)
        for idx, cov in node_inst.to_edges.iteritems():
            to_node = self.G[idx]
            to_node.from_edges.pop(node_inst.idx, None)
    
    def removeEdge(self, fro_idx, to_idx):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ("Node indexes should be INT", fro_idx, to_idx)
            return False
        if fro_idx in self.G and to_idx in self.G:            
            self.G[fro_idx].to_edges.pop(to_idx, None) # Remove edge from previous node            
            self.G[to_idx].to_edges.pop(fro_idx, None) # Remove edge to next node
        else:
            print ("removeEdge(): Nodes not exist", fro_idx, to_idx)
            return False
        return True
    
    def lookupMutationGraphNode(self, original_seq_len, node_from, node_to, snp_char, cov, node_type):
        if self.graph_node_len == original_seq_len:
            # Graph does not have any mutation node so we add
            new_ind = self.graph_node_len
            self.addNode(new_ind, snp_char, cov, node_type)
            self.addEdge(node_from, new_ind, cov)
            self.addEdge(new_ind, node_to, cov)
            return True
        elif self.graph_node_len < original_seq_len:
            print ('graph structure was not valid upon initiation')
            return False
        for idx, node in self.G.iteritems():
            if idx < original_seq_len:
                continue
            if node_from in node.from_edges and node_to in node.to_edges and node.base == snp_char:
                node.cov += cov
                self.addEdge(node_from, idx, cov)
                self.addEdge(idx, node_to, cov)
                return True
        # We don't have identical existing mutation node so we add
        new_ind = self.graph_node_len
        self.addNode(new_ind, snp_char, cov, node_type)
        self.addEdge(node_from, new_ind, cov)
        self.addEdge(new_ind, node_to, cov)
        return True
    
    def refNodeExists(self, ref_coor, ref_base, prev_ind, next_ind):
        if prev_ind not in self.G or next_ind not in self.G:
            return False
        for i in range(0,len(ref_base)):
            if ref_coor+i not in self.G:
                return False
            if self.G[ref_coor+i].base != ref_base[i]:
                print ('VCF coordinate does not match with reference',\
                ref_coor, self.G[ref_coor+i].base, ref_base)
                return False
        return True
        
    def addSNP(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor - 1
        next_ind = ref_coor + 1
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base, cov, 'SNP', prev_ind, next_ind, ref_coor, cov)
        return True

    def addINS(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor
        next_ind = ref_coor + 1
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base[1:], cov, 'INS', prev_ind, next_ind, ref_coor, cov)
        return True

    def addSUB(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor - 1
        next_ind = ref_coor + len(mut_base)
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base, cov, 'SUB', prev_ind, next_ind, ref_coor, cov)
        return True

    def addDEL(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor
        next_ind = ref_coor + len(ref_base)
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation('', cov, 'DEL', prev_ind, next_ind, ref_coor, cov)
        return True    
    