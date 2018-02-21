# -*- coding: utf-8 -*-
"""
Created on Tue Aug 01 20:13:22 2017

@author: swoo
"""
import parameter


class genTestFASTQ(object):
    def __init__(self):
        self.prmtr = parameter.Parameter()
        self.gff = 'F:\\Ensembl\\out_seq.txt'
        self.out_FASTQ = 'F:\\Ensembl\\out_seq.fq'
        self.line1_str = '@ZZZ'
        self.line3_str = '+ZZZ\n'
        self.line5_str = '???\n'
        self.cnt = 0
        self.Process()
        
    def Process(self):
        self.parseGFF()

    def parseGFF(self):
        trans_seq = ''
        list_coor_tree = []
        with open(self.gff,'r') as in_gff, open(self.out_FASTQ,'w') as oFASTQ:
            for line in in_gff:                
                if line.startswith('A') or line.startswith('T') or line.startswith('C') or line.startswith('G'):
                    trans_seq = line.strip()
                else:
                    for i in range(len(trans_seq)-self.prmtr.max_idx):
                        self.writeSeq(trans_seq[i:i+self.prmtr.max_idx], oFASTQ)
                    self.cnt += 1
                        
    def writeSeq(self, seq, oFile_inst):
        oFile_inst.write(self.line1_str + str(self.cnt) + '\n')
        oFile_inst.write(seq + '\n')
        oFile_inst.write(self.line3_str)
        oFile_inst.write(self.line5_str)
                        
if __name__=='__main__':
    obj = genTestFASTQ()
    
                        
                        
    