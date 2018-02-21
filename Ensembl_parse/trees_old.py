# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 20:41:08 2017

@author: Sunghee Woo
"""
from collections import defaultdict
import pickle

class Node:
      def __init__(self, base_str, coverage, node_type):
           #self.idx = idx
           self.base = base_str
           self.cov = coverage
           self.from_edges = defaultdict(int)
           self.to_edges = defaultdict(int)
           self.visited = 0
           self.type = node_type
           
class TGraph:
    def __init__(self):
        self.G = {}
        self.graph_node_len = 0
        
    def Process(self):
        print ('max=')
        
    def addNode(self, idx, base_str, coverage, node_type):
        try:
            idx = int(idx)
        except ValueError:
            #print "Node idx should be INT", idx
            return False
        if idx in self.G:
            #print "Node already exist, adding coverage to existing node", idx, base_str, coverage
            self.G[idx].cov += coverage
        else:
            self.G[idx] = Node(base_str, coverage, node_type)
            self.graph_node_len += 1
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
            print ("Node indexes should be INT", fro_idx, to_idx)
            return False
        if fro_idx in self.G and to_idx in self.G:            
            self.G[fro_idx].to_edges.pop(to_idx, None) # Remove edge from previous node            
            self.G[to_idx].to_edges.pop(fro_idx, None) # Remove edge to next node
        else:
            print ("removeEdge(): Nodes not exist", fro_idx, to_idx)
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

class Parse(object):
    def __init__(self, gff, fasta):
        self.genes = {}
        self.gff = gff
        self.fasta = fasta
        self.chr = ''
        self.seq = ''
        self.parseDNA()
        print ('parsed DNA')
        self.parseGFF()
        
    def parseDNA(self):
        with open(self.fasta,'r') as in_fasta:#, open("F:\\Ensembl\\dna_chr1.txt",'w') as o_fasta:
            header = in_fasta.readline()
            self.chr = header.split(' ')[0]
            seqlist = in_fasta.read().splitlines()
            self.seq = ''.join(seqlist)
            #pickle.dump(self.seq, open("F:\\Ensembl\\dna_chr1.p", "wb"))
            #o_fasta.write(self.seq)
#            for line in in_fasta:
#                if line.startswith('>'):
#                    self.chr = line[1:].split(' ')[0]
#                else:
#                    self.seq += line.strip()
            
    def parseGFF(self):
        in_transcript = False
        trans_id_cds = -1
        with open(self.gff,'r') as in_gff:
            for line in in_gff:
                if line.startswith('#'):
                    continue
                data = line.split('\t')
                if data[0] != '1':
                    continue
                #if data[2] == 'transcript':
                if data[2] == 'mRNA':
#                    if trans_id_cds != -1:
#                        print (trans_id_cds, trans_seq)
                    trans_seq = ''
                    #chr_id = data[0]
                    tags = data[8].split(';')
                    for i in tags:
                        if i.startswith('Parent=gene:'):
                            gene_id = i.replace('Parent=gene:','')
                        if i.startswith('ID=transcript:'):
                            trans_id = i.replace('ID=transcript:','')
                    if gene_id not in self.genes:
                        self.genes[gene_id] = TGraph()
                    in_transcript = True
                    node_from = -1
                    continue
                elif data[2] == 'exon':
                    continue
                elif data[2].find('UTR') > -1:
                    continue
                elif data[2] == 'CDS' and in_transcript == True:
                    tags = data[8].split(';')
                    for i in tags:
#                        if i.startswith('Parent=gene:'):
#                            gene_id = i.replace('Parent=gene:','')
                        if i.startswith('Parent=transcript:'):
                            trans_id_cds = i.replace('Parent=transcript:','')
                            break
                    if trans_id != trans_id_cds:
                        print ('gff3 parsing is going wrong', trans_id, trans_id_cds)                    
                    start = int(data[3])-1
                    end = int(data[4])
                    seq = self.getSeq(start, end)
                    for i in range(0,len(seq)):
                        if self.genes[gene_id].G == {}:
                            self.genes[gene_id].addNode(-1, 'S', 0, 'source')
                            node_from = -1
                        node_idx = start + i
                        self.genes[gene_id].addNode(node_idx, seq[i], 0, trans_id_cds)
                        self.genes[gene_id].addEdge(node_from, node_idx, 0)
                        node_from = node_idx
                    trans_seq += seq
                else:
                    in_transcript = False
                    continue
        pickle.dump(self.genes, open("F:\\Ensembl\\genes_chr1.p", "wb"))
        
        
    def getSeq(self, start, end):
        return self.seq[start:end]
    

class SuffixTree(object):
    
    class Node(object):
        def __init__(self, lab):
            self.lab = lab # label on path leading to this node
            self.out = {}  # outgoing edges; maps characters to nodes
    
    def __init__(self, s):
        """ Make suffix tree, without suffix links, from s in quadratic time
            and linear space """
        s += '$'
        self.root = self.Node(None)
        self.root.out[s[0]] = self.Node(s) # trie for just longest suf
        # add the rest of the suffixes, from longest to shortest
        for i in xrange(1, len(s)):
            # start at root; we’ll walk down as far as we can go
            cur = self.root
            j = i
            while j < len(s):
                if s[j] in cur.out:
                    child = cur.out[s[j]]
                    lab = child.lab
                    # Walk along edge until we exhaust edge label or
                    # until we mismatch
                    k = j+1 
                    while k-j < len(lab) and s[k] == lab[k-j]:
                        k += 1
                    if k-j == len(lab):
                        cur = child # we exhausted the edge
                        j = k
                    else:
                        # we fell off in middle of edge
                        cExist, cNew = lab[k-j], s[k]
                        # create “mid”: new node bisecting edge
                        mid = self.Node(lab[:k-j])
                        mid.out[cNew] = self.Node(s[k:])
                        # original child becomes mid’s child
                        mid.out[cExist] = child
                        # original child’s label is curtailed
                        child.lab = lab[k-j:]
                        # mid becomes new child of original parent
                        cur.out[s[j]] = mid
                else:
                    # Fell off tree at a node: make new edge hanging off it
                    cur.out[s[j]] = self.Node(s[j:])
    
    def findPath(self, s):
        """ Follow path given by s.  If we fall off tree, return None.  If we
            finish mid-edge, return (node, offset) where 'node' is child and
            'offset' is label offset.  If we finish on a node, return (node,
            None). """
        cur = self.root
        i = 0
        while i < len(s):
            c = s[i]
            if c not in cur.out:
                return (None, None) # fell off at a node
            child = cur.out[s[i]]
            lab = child.lab
            j = i+1
            while j-i < len(lab) and j < len(s) and s[j] == lab[j-i]:
                j += 1
            if j-i == len(lab):
                cur = child # exhausted edge
                i = j
            elif j == len(s):
                return (child, j-i) # exhausted query string in middle of edge
            else:
                return (None, None) # fell off in the middle of the edge
        return (cur, None) # exhausted query string at internal node
    
    def hasSubstring(self, s):
        """ Return true iff s appears as a substring """
        node, off = self.findPath(s)
        return node is not None
    
    def hasSuffix(self, s):
        """ Return true iff s is a suffix """
        node, off = self.findPath(s)
        if node is None:
            return False # fell off the tree
        if off is None:
            # finished on top of a node
            return '$' in node.out
        else:
            # finished at offset 'off' within an edge leading to 'node'
            return node.lab[off] == '$'
            
        
if __name__=='__main__':
    path = Parse('F:\\Ensembl\\Homo_sapiens.GRCh38.88.chr.gff3','F:\\Ensembl\\Homo_sapiens.GRCh38.dna_rm.chromosome.1.fa')
        
        
        
        