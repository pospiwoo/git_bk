# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 20:56:29 2017

@author: Sunghee Woo
"""
class Parameter(object):
    def __init__(self):
        self.trie_flag = True # True False
        self.trie_flag = False # True False
        self.max_idx = 110
        self.min_idx = 9
        self.max_mismatch = 0
        self.gff = 'F:\\Ensembl\\out_seq.txt'
        self.gff = 'F:\\Ensembl\\Homo_sapiens.GRCh38.88.chr_6.gff3' # 'F:\\Ensembl\\Homo_sapiens.GRCh38.dna_rm.chromosome.1.fa'
#        self.gff = 'F:\\Ensembl\\out_seq_2.txt'
        self.fasta = 'F:\\Ensembl\\unmasked\\Homo_sapiens.GRCh38.dna.chromosome.1.fa'
        self.vcf = 'F:\\Ensembl\\common_all_20170403_6.vcf'
        self.fastq = 'F:\\sra\\ERR904420.fastq'
        self.fastq = 'F:\\Ensembl\\out_seq.fq'
        self.gen_rand_error = True # False True
        self.debug_trans_seq = '' # 'F:\\Ensembl\\out_seq.txt'
        
        