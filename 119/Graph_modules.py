# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 12:47:31 2016

@author: Sunghee Woo
"""
import os, math #, re, string, sys
from collections import defaultdict
import matplotlib.pyplot as plt
import Parameters as param

class Node:
      def __init__(self, idx, base_str, coverage, node_type):
           self.idx = idx
           self.base = base_str
           self.cov = coverage
           self.from_edges = defaultdict(float)
           self.to_edges = defaultdict(float)
           self.visited = 0
           self.type = node_type
           
class SpliceGraph:
    def __init__(self, G_view_inst):
        self.G = {}
        self.GraphView = G_view_inst
        self.graph_node_len = 0
        self.max_seq = []
        self.max_seq_types = []
        self.max_seq_cov = []
        self.max_seq_score = []
        self.DynamicP_track = []
        self.DynamicP_score = []
        self.ideal_barcode_cov = []
        
    def Process(self):
        maxa = self.findLongestPath('a')
        print 'max=', maxa
        
    def addNode(self, idx, base_str, coverage, node_type):
        try:
            idx = int(idx)
        except ValueError:
            print "Node idx should be INT", idx
            return False
        if idx in self.G:
            print "Node already exist, adding coverage to existing node", idx, base_str, coverage
            self.G[idx].cov += coverage
        else:
            self.G[idx] = Node(idx, base_str, coverage, node_type)
            self.graph_node_len += 1
        return True
        
    def addEdge(self, fro_idx, to_idx, cov):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print "Node indexes should be INT", fro_idx, to_idx
            return False
        if fro_idx in self.G and to_idx in self.G:            
            self.G[fro_idx].to_edges[to_idx] += cov # Assign edge from previous node            
            self.G[to_idx].from_edges[fro_idx] += cov # Assign edge to next node
        else:
            print "Nodes not exist", fro_idx, to_idx
            return False
        return True
        
    def removeNode(self, node_inst):
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
            print "Node indexes should be INT", fro_idx, to_idx
            return False
        if fro_idx in self.G and to_idx in self.G:            
            self.G[fro_idx].to_edges.pop(to_idx, None) # Remove edge from previous node            
            self.G[to_idx].to_edges.pop(fro_idx, None) # Remove edge to next node
        else:
            print "removeEdge(): Nodes not exist", fro_idx, to_idx
            return False
        return True
        
    def lookupMutationNode(self, original_seq_len, node_from, node_to, snp_char, cov, node_type):
        if self.graph_node_len == original_seq_len:
            # Graph does not have any mutation node so we add
            new_ind = self.graph_node_len
            self.addNode(new_ind, snp_char, cov, node_type)
            self.addEdge(node_from, new_ind, cov)
            self.addEdge(new_ind, node_to, cov)
            return True
        elif self.graph_node_len < original_seq_len:
            print 'graph structure was not valid upon initiation'
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
    
    def initDynamicPath(self):
        self.DynamicP_track = []
        self.DynamicP_score = []
        for idx, node in self.G.iteritems():
            self.DynamicP_track.append(-1)
            self.DynamicP_score.append(0)
        
    def findDynamicPath(self, curr_idx):
        curr_node = self.G[curr_idx]
        curr_node.visited = 1
        if not curr_node.to_edges:
            self.DynamicP_track[curr_node.idx] = -2
            self.DynamicP_score[curr_node.idx] = curr_node.cov
            return
        for to_idx in curr_node.to_edges:
            next_node = self.G[to_idx]
            if next_node.visited == 0:
                self.findDynamicPath(to_idx)
                new_dist = next_node.cov + curr_node.cov                
                if new_dist > self.DynamicP_score[curr_node.idx]:
                    self.DynamicP_track[curr_node.idx] = next_node.idx
                    self.DynamicP_score[curr_node.idx] = new_dist
                    #print "curr_idx", curr_idx, to_idx
        curr_node.visited = 0
        return
        
    def findBruteForcePath(self, curr_idx):
        sum_dist = 0
        max_dist = 0
        curr_node = self.G[curr_idx]
        curr_node.visited = 1
        for to_idx in curr_node.to_edges:
            next_node = self.G[to_idx]
            if next_node.visited == 0:
                sum_dist = next_node.cov + self.findDynamicPath(to_idx)
                if sum_dist > max_dist:
                    max_dist = sum_dist
                    #print "curr_idx", curr_idx, to_idx
                    self.DynamicP_track[curr_idx] = to_idx
                    self.DynamicP_score[curr_idx] = max_dist + curr_node.cov
        curr_node.visited = 0
        return max_dist
    
    def greedyPath(self, curr_idx):
        curr_node = self.G[curr_idx]
        curr_node.visited = 1
        max_cov = -1.0
        next_idx = -1
        self.max_seq.append(curr_node.base)
        self.max_seq_types.append(self.GraphView.checkMutationString(curr_node))                
        if len(curr_node.to_edges) > 0:
            for to_idx in curr_node.to_edges:
                next_node = self.G[to_idx]
                if next_node.visited == 0 and next_node.cov > max_cov:
                    max_cov = next_node.cov
                    next_idx = to_idx
            self.greedyPath(next_idx)
        return
    
    def greedyQualityPath(self, curr_idx):
        curr_node = self.G[curr_idx]
        curr_node.visited = 1
        max_cov = -1.0
        next_idx = -1
        self.max_seq_types.append(self.GraphView.checkMutationString(curr_node))
        self.max_seq.append(curr_node.base)
        #self.max_seq_score.append(str(int(min(math.floor(curr_node.cov+0.5),9.0)))) # tmp score
        self.max_seq_score.append(str(int(min(math.floor(curr_node.cov),9.0)))) # tmp score
        self.max_seq_cov.append(str(curr_node.cov))
        if len(curr_node.to_edges) > 0:
            for to_idx in curr_node.to_edges:
                next_node = self.G[to_idx]
                if next_node.visited == 0 and next_node.cov > max_cov:
                    max_cov = next_node.cov
                    next_idx = to_idx
            self.greedyQualityPath(next_idx)
        return
        
    def GraphTrimming(self, threshold):
        #3.315326811 SNP, 2.401673047 SNP, 1.865305544 SNP
        for idx, node in self.G.iteritems():
            #print node.idx, node.base, node.from_edges, node.to_edges, node.cov
            if node.type != 'ref' and node.cov < float(threshold):
                self.removeNode(node)
                
    def CallInsertions(self, threshold):
        for idx, node in self.G.iteritems():
            if node.type == 'INS':
                if len(node.from_edges) != 1 or len(node.to_edges) != 1:
                    print "CallInsertions(): Insertion node must have only one fromEdge and toEedge", len(node.from_edges), node.from_edges, node.to_edges
                    return False
                for from_edge in node.from_edges:
                    from_e = from_edge
                for to_edge in node.to_edges:
                    to_e = to_edge                    
                if self.G[from_e].type != 'ref' or self.G[to_e].type != 'ref':
                    print "CallInsertions(): Insertion node must have only reference edges connected"
                    return False
                avg_threshold = (self.G[from_e].cov + self.G[to_e].cov) / 2.0 / 2.0 # avg
                if node.cov >= avg_threshold and node.cov >= float(threshold):
                    self.removeEdge(from_e, to_e)


class GraphView:
    def __init__(self, path_inst):
        self.path = path_inst
            
    def print_matrix(self, matrix):
        for mat_line in matrix:
            print('\t'.join(map(str,mat_line)))
        return
        
    def PrintGraph(self, Graph_inst, out_file_name, origin_gene):
        out_file_name = os.path.join(self.path.out_dir, out_file_name)
        with open(out_file_name,'w') as oFile: # print graph structure
            oFile.write('>'+origin_gene+'\n')
            oFile.write(''.join(Graph_inst.max_seq)+'\n')
            oFile.write(''.join(Graph_inst.max_seq_types)+'\n')
            for idx, node in Graph_inst.G.iteritems():
                #print node.idx, node.base, node.from_edges, node.to_edges, node.cov
                oFile.write(str(node.idx) + '\t' +
                    str(node.base) + '\t' +
                    'from:' +
                    str(node.from_edges) + '\t' +
                    'to:' +
                    str(node.to_edges) + '\t' +
                    str(node.type) + '\t' +
                    str(node.cov) + '\n' )
        
    def PrintGraphAvg(self, Graph_inst, out_file_name, div, origin_gene):
        out_file_name = os.path.join(self.path.out_dir, out_file_name)        
        newList = [float(x) / float(div) for x in Graph_inst.max_seq_cov]
        plt.plot(newList)
        plt.ylabel('average coverage')
        plt.savefig(out_file_name.replace('.txt','.pdf'))
        plt.close()
        with open(out_file_name,'w') as oFile: # print graph structure
            oFile.write('>'+origin_gene+'\n')
            oFile.write(''.join(Graph_inst.max_seq)+'\n')
            oFile.write(''.join(Graph_inst.max_seq_types)+'\n')
            for idx, node in Graph_inst.G.iteritems():
                node = Graph_inst.G[idx]
                #print node.idx, node.base, node.from_edges, node.to_edges, node.cov
                oFile.write(str(node.idx) + '\t' +
                    str(node.base) + '\t' +
                    'from:' +
                    str(node.from_edges) + '\t' +
                    'to:' +
                    str(node.to_edges) + '\t' +
                    str(node.type) + '\t' +
                    str(node.cov/div) + '\n' )
        
    def PrintGraphAvgWithIdeal(self, Graph_inst, out_file_name, div, ideal):
        out_file_name = os.path.join(self.path.out_dir, out_file_name)
        newList = [float(x) / float(div) for x in Graph_inst.max_seq_cov]
        plt.plot(ideal, label="perfect cov")
        plt.plot(newList, label="exp cov")
        plt.ylabel('average coverage')
        plt.xlabel('base #')
        #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
        plt.legend()
        plt.savefig(out_file_name.replace('.txt','.pdf'))
        plt.close()
        with open(out_file_name,'w') as oFile: # print graph structure
            oFile.write('>\n')
            oFile.write(''.join(Graph_inst.max_seq)+'\n')
            oFile.write(''.join(Graph_inst.max_seq_types)+'\n')
            for idx, node in Graph_inst.G.iteritems():
                node = Graph_inst.G[idx]
                #print node.idx, node.base, node.from_edges, node.to_edges, node.cov
                oFile.write(str(node.idx) + '\t' +
                    str(node.base) + '\t' +
                    'from:' +
                    str(node.from_edges) + '\t' +
                    'to:' +
                    str(node.to_edges) + '\t' +
                    str(node.type) + '\t' +
                    str(node.cov/div) + '\n' )
        
    def PrintGraphMaxPath(self, Graph_inst, out_file_name, origin_gene):
        out_file_name = os.path.join(self.path.out_dir, out_file_name)
        with open(out_file_name,'w') as oFile: # print graph structure
            oFile.write('>'+origin_gene+'\n')
            for idx in range(len(Graph_inst.max_seq)):
                if param.mask_global_reads() and param.mask_global_N_thres() > float(Graph_inst.max_seq_cov[idx]):
                    oFile.write('N')
                else:
                    oFile.write(Graph_inst.max_seq[idx])
            oFile.write('\n')                
            for idx in range(len(Graph_inst.max_seq)):
                if param.mask_global_reads() and param.mask_global_N_thres() > float(Graph_inst.max_seq_cov[idx]):
                    oFile.write('u')
                else:
                    oFile.write(Graph_inst.max_seq_types[idx])
            oFile.write('\n')            
            for idx in range(len(Graph_inst.max_seq)):
                oFile.write(str(idx) + '\t')
                oFile.write(Graph_inst.max_seq[idx] + '\t')
                oFile.write(Graph_inst.max_seq_types[idx] + '\t')
                oFile.write(Graph_inst.max_seq_cov[idx] + '\t')
                oFile.write(Graph_inst.max_seq_score[idx] + '\n')
                
    def PrintFullReads(self, Graph_inst, oFile, header_str, gene_str):
        if header_str:
            oFile.write('>'+gene_str+';'+header_str+'\n')
        else:
            oFile.write('>'+gene_str+'\n')
        oFile.write(''.join(Graph_inst.max_seq)+'\n')
        oFile.write(''.join(Graph_inst.max_seq_types)+'\n')
        if Graph_inst.max_seq_score != []:
            oFile.write(''.join(Graph_inst.max_seq_score)+'\n')
        
    def PrintFullReadsDynamicP(self, Graph_inst, oFile, header_str, gene_str):
        if header_str:
            oFile.write('>'+header_str+':'+gene_str+'\n')
        else:
            oFile.write('>'+gene_str+'\n')
        full_read = ''
        full_type = ''
        full_score = []
        next_idx = 0
        while next_idx != -2:
            curr_node = Graph_inst.G[next_idx]
            full_type += self.checkMutationString(curr_node)
            full_read += curr_node.base
            full_score.append(str(int(min(math.floor(curr_node.cov+0.5),9.0)))) # tmp score
            next_idx = Graph_inst.DynamicP_track[next_idx]
        oFile.write(full_read + '\n')
        oFile.write(full_type + '\n')
        oFile.write(''.join(full_score) + '\n')
#        oFile.write(','.join(map(str,Graph_inst.DynamicP_track))+'\n')
#        oFile.write(','.join(map(str,Graph_inst.DynamicP_score))+'\n')

    def checkMutationString(self, node):
        return_str = ''
        if float(node.cov) < param.mask_ref_N_thres() and node.type == 'ref':
            if param.mask_ref_reads():
                node.base = 'N'
            return 'u'
        if float(node.cov) < param.mask_snp_N_thres() and node.type == 'SNP':
            if param.mask_snp_reads():
                node.base = 'N'
            return 'u'
        if float(node.cov) < param.mask_mu_N_thres() and node.type != 'ref':
            if param.mask_mu_reads():
                node.base = 'N'
            return 'u'
        if node.type == 'ref':
            return_str = '-'
        elif node.type == 'SNP':
            return_str = '*'
        elif node.type == 'SUB':
            return_str = '^'
        elif node.type == 'INS':
            return_str = '+'
        elif node.type == 'DEL':
            return_str = '='
        else:
            return_str = '?'
        return return_str
        

#if __name__=='__main__':
#    Graph_view = GraphView()
#    Graph = SpliceGraph(Graph_view)
#    Graph.Process()

