# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 20:04:58 2017

@author: Sunghee Woo
"""
import os, sys, pickle
import graph, trie, tree, parameter
#sys.setrecursionlimit(30000)

class ConstructGraph(object):
    def __init__(self):
        self.prmtr = parameter.Parameter()
        self.GGraph = graph.GGraph()
        self.GGraph.addNode(-1, 's', 0, 'source')
        self.chr = ''
        self.seq = ''
        self.Process()
        
    def Process(self):
        print ('Reading reference DNA')
        self.parseDNA()
        print ('Creating transcript graph from GFF')
        self.parseGFF()
        print ('Creating variant graph fvrom VCF')
        self.parseVCF()
        
    def parseDNA(self):
        with open(self.prmtr.fasta,'r') as in_fasta:#, open("F:\\Ensembl\\dna_chr1.txt",'w') as o_fasta:
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
        tmp_cnt = 0
        if self.prmtr.trie_flag:
            self.tree = trie.SuffixTrie(self.prmtr)
        else:
            self.tree = tree.SuffixTree(self.prmtr)
        trans_seq = ''
        list_coor_tree = []
        with open(self.prmtr.gff,'r') as in_gff:
            if self.prmtr.debug_trans_seq != '':
                debug_seq = open('F:\\Ensembl\\out_seq.txt','w')
            for line in in_gff:
                if line.startswith('#'):
                    continue
                data = line.split('\t')
                if data[0] != '1':
                    continue
                if data[2] == 'mRNA':
                    #chr_id = data[0]
                    tags = data[8].split(';')
                    for i in tags:
                        if i.startswith('Parent=gene:'):
                            gene_id = i.replace('Parent=gene:','')
                        if i.startswith('ID=transcript:'):
                            trans_id = i.replace('ID=transcript:','')
                    in_transcript = True
                    node_from = -1
                    tmp_cnt += 1
                    if list_coor_tree != []:
                        self.tree.makeST(self.tree, trans_seq, list_coor_tree)
                        if self.prmtr.debug_trans_seq != '':
                            debug_seq.write(trans_seq+'\n'+','.join(str(x) for x in list_coor_tree)+'\n')
                        list_coor_tree = []
                        trans_seq = ''
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
                        node_idx = start + i
                        self.GGraph.addNode(node_idx, seq[i], 0, trans_id_cds)
                        self.GGraph.addEdge(node_from, node_idx, 0)
                        node_from = node_idx
                        list_coor_tree.append(node_idx)
                    trans_seq += seq
                else:
                    in_transcript = False
                    continue
#        pickle.dump(self.GGraph, open("F:\\Ensembl\\genes_chr1.p", "wb"))
        pickle.dump(self.tree, open("F:\\Ensembl\\trie_chr1_6.p", "wb"))
        if self.prmtr.debug_trans_seq != '':
            debug_seq.close()
        print ('num genes:', tmp_cnt)
        print ('G len ref:', self.GGraph.graph_node_len)
        
    def parseVCF(self):
        with open(self.prmtr.vcf,'r') as in_vcf:
            for line in in_vcf:
                if line.startswith('#'):
                    continue
                data = line.strip().split('\t')
                if data[0] != '1':
                    continue
                chr_id = data[0]
                start = int(data[1])-1
                ref = data[3]
                mut_list = data[4].split(',')
                for mut in mut_list:
                    if len(ref) == 1: # SNP / insertion
                        if len(mut) == 1: # SNP
                              self.GGraph.addSNP(start, ref, mut, 0)
                        else: # insertion
                              self.GGraph.addINS(start, ref, mut, 0)
                    elif len(ref) > 1: # deletion / substitution
                        if len(mut) == len(ref): # substitution
                              self.GGraph.addSUB(start, ref, mut, 0)
                        else: # deletion
                              self.GGraph.addDEL(start, ref, mut, 0)
                    else:
                        print ('unexpected format in VCF:', line)
                        
        pickle.dump(self.GGraph, open("F:\\Ensembl\\genes_mu_chr1_6.p", "wb"))
        print ('G len:', self.GGraph.graph_node_len)
        print ('SNP:', len(self.GGraph.SNP))
        print ('SUB:', len(self.GGraph.SUB))
        print ('INS:', len(self.GGraph.INS))
        print ('DEL:', len(self.GGraph.DEL))        
        
    def getSeq(self, start, end):
        return self.seq[start:end]    


class unpickleGraph(object):
    def __init__(self):
        self.prmtr = parameter.Parameter()
        print ('loading gene pickle')
        self.GGraph = pickle.load( open("F:\\Ensembl\\genes_mu_chr1_6.p", "rb" ) )
        print ('loading tree pickle')
        self.tree = pickle.load( open("F:\\Ensembl\\trie_chr1_6.p", "rb" ) )
        print ('pickle loaded')


class testConstructGraph(object):
    def __init__(self):
        self.prmtr = parameter.Parameter()
        self.GGraph = graph.GGraph()
        self.GGraph.addNode(-1, 's', 0, 'source')
        self.Process()
        
    def Process(self):
        self.parseGFF()

    def parseGFF(self):
        cnt = 0
        if self.prmtr.trie_flag:
            self.tree = trie.SuffixTrie(self.prmtr)
        else:
            self.tree = tree.SuffixTree(self.prmtr)
        trans_seq = ''
        list_coor_tree = []
        with open(self.prmtr.gff,'r') as in_gff:
            for line in in_gff:                
                if line.startswith('A') or line.startswith('T') or line.startswith('C') or line.startswith('G'):
                    trans_seq = line.strip()
                else:
                    list_coor_tree = [int(x) for x in line.strip().split(',')]
                    self.tree.makeST(self.tree, trans_seq, list_coor_tree[:len(trans_seq)])
                    cnt += 1
                    if cnt % 50 == 0:
                        print (cnt, 'trans indexed')
    
if __name__=='__main__':
    parse_obj = ConstructGraph()
#    parse_obj.Process()
        