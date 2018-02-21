# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 09:02:03 2016

@author: Sunghee Woo
"""
#import os, sys, re, operator, numarray, imap, unittest
class Alignment():
    def __init__(self):
        self.match = 2 #1
        self.mismatch = -1 #-1
        self.gap = -30 #-1
        self.target_seq = ''
        self.sixmer = ''

    def ProcessHamming(self, molecule_inst, target_seq, sixmer):
        self.target_seq = target_seq
        self.sixmer = sixmer
        snp_dic = {}
        for i in xrange(0,len(target_seq)-6):
            target_chunk = target_seq[i:i+6]
            diffs, diff_idx, diff_char = self.hammingDist(target_chunk, sixmer)
            if diffs == 1:
                snp_dic[diff_idx] = [diff_char, target_chunk]
        return snp_dic

    def hammingDist(self, target_str, six_str):
        diffs = 0
        if len(target_str) == len(six_str):
            for i in xrange(len(target_str)):
                if target_str[i] != six_str[i]:
                    diffs += 1
                    last_diff_idx = i
                    diff_char = six_str[i]
        else:
            print "hammingDist(): length of strings are different", target_str, six_str
        return diffs, last_diff_idx, diff_char

    def ProcessSW(self, molecule_inst, target_seq, sixmer):
        self.target_seq = target_seq
        self.sixmer = sixmer
        self.num_rows = len(self.target_seq) + 1
        self.num_cols = len(self.sixmer) + 1
        self.scoring_table = [[0 for col in range(self.num_cols)] for row in range(self.num_rows)]
        self.start_pos   = None
        
        # Scoring table
        self.create_scoring_table()
        
        # Table DP_table_traceback
        self.target_seq_aligned, self.sixmer_aligned = self.DP_table_traceback(molecule_inst)
        #self.target_seq_aligned, self.sixmer_aligned = self.DP_table_traceback_diag()
        
        assert len(self.target_seq_aligned) == len(self.sixmer_aligned), 'aligned strings are not the same size'

        #for i in xrange(0,len(self.scoring_table)):
        #    print ", ".join(str(self.scoring_table[i])  )
        #print self.target_seq_aligned, self.sixmer_aligned, '6mer', self.sixmer
        if not self.target_seq_aligned.find('-') > -1 and \
            not self.sixmer_aligned.find('-') > -1 and \
            len(self.target_seq_aligned) == len(self.sixmer_aligned) and \
            len(self.target_seq_aligned) == 6:
            #for i in xrange(0,len(self.target_seq_aligned)):
            #    if self.target_seq_aligned[i] == self.sixmer_aligned[i]:
            #print self.target_seq_aligned, self.sixmer_aligned, '6mer', self.sixmer
            return self.target_seq_aligned, self.sixmer_aligned
        else:
            return '', ''      
            
    def create_scoring_table(self):
        max_score = 0
        max_pos   = None
        for i in range(1, self.num_rows):
            for j in range(1, self.num_cols):
                score = self.score_lookup(i, j)
                if score > max_score:
                    max_score = score
                    max_pos   = (i, j)
        
                self.scoring_table[i][j] = score
        assert max_pos is not None, 'the x, y position with the highest score was not found'
        self.start_pos = max_pos
        
    def score_lookup(self, x, y):
        similarity = self.match if self.target_seq[x - 1] == self.sixmer[y - 1] else self.mismatch
        diag_score = self.scoring_table[x - 1][y - 1] + similarity
        up_score   = self.scoring_table[x - 1][y] + self.gap
        left_score = self.scoring_table[x][y - 1] + self.gap
        return max(0, diag_score, up_score, left_score)

    def DP_table_traceback(self, molecule):
        # DIAG: Match/Mismatch
        # UP: Gap in target
        # LEFT: Gap in sixmer
        END, DIAG, UP, LEFT = range(4)
        self.aligned_target_seq = []
        self.aligned_sixmer = []
        x, y         = self.start_pos
        move         = self.next_move(x, y)
        while move != END:
            if move == DIAG:
                self.aligned_target_seq.append(self.target_seq[x - 1])
                self.aligned_sixmer.append(self.sixmer[y - 1])
                #if self.target_seq[x - 1] != self.sixmer[x - 1]:
                #    molecule.
                x -= 1
                y -= 1
            elif move == UP:
                self.aligned_target_seq.append(self.target_seq[x - 1])
                self.aligned_sixmer.append('-')
                x -= 1
            else:
                self.aligned_target_seq.append('-')
                self.aligned_sixmer.append(self.sixmer[y - 1])
                y -= 1

            move = self.next_move(x, y)

        self.aligned_target_seq.append(self.target_seq[x - 1])
        self.aligned_sixmer.append(self.sixmer[y - 1])
        
        return ''.join(reversed(self.aligned_target_seq)), ''.join(reversed(self.aligned_sixmer))

    def next_move(self, x, y):
        diag = self.scoring_table[x - 1][y - 1]
        up   = self.scoring_table[x - 1][y]
        left = self.scoring_table[x][y - 1]
        if diag >= up and diag >= left:     # Tie-diag
            return 1 if diag != 0 else 0    # 1-diag, 0-end
        elif up > diag and up >= left:      # Tie-up
            return 2 if up != 0 else 0      # up or end
        elif left > diag and left > up:
            return 3 if left != 0 else 0    # left or end
        else:
            raise ValueError('invalid move during DP_table_traceback')

    def print_matrix(matrix):
        for row in matrix:
        	for col in row:
        		print('{0:>4}'.format(col))
        	print()
         

if __name__=='__main__':
	target_seq = 'QTTATCCTTGTCCTTGTCCTA'
	sixmer = "GGATTC"	
	Alignment = Alignment()
	#Alignment.Process(target_seq, sixmer)
