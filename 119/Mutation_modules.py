# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:37:47 2016

@author: Sunghee Woo
"""
import os
import Graph_modules as Graph
import Parameters as param

class MutationModel:
    def __init__(self, Assm_model, Graph_view):
        self.default_cov = param.default_mutation_cov()
        self.default_ref_cov = param.default_ref_cov()
        self.Assm_model = Assm_model
        self.Graph_view = Graph_view
        self.mutation_file_name = Assm_model.mutation_file_name
        self.sixmer_dic_mutation = {}
        self.sixmer_dic_l_mutation = {}
        self.Graph = {}
        self.mutations = {}
        self.initGraph()
        self.NN_subst = ['AA', 'AT', 'AC', 'AG',\
                         'TA', 'TT', 'TC', 'TG',\
                         'CA', 'CT', 'CC', 'CG',\
                         'GA', 'GT', 'GC', 'GG']
        
    def initGraph(self):
        for gene_str, seq in self.Assm_model.target_sequences.iteritems():
            self.Graph[gene_str] = Graph.SpliceGraph(self.Graph_view)
            for i in xrange(len(seq)):
                self.Graph[gene_str].addNode(i, seq[i], self.default_ref_cov, 'ref')
            for i in xrange(1, len(seq)-1):
                self.Graph[gene_str].addEdge(i-1, i, self.default_ref_cov)
                self.Graph[gene_str].addEdge(i, i+1, self.default_ref_cov)
        if param.enable_known_mu():
            self.parseMutationFile()
            self.addMutations()
        self.Assm_model.Assm_view.print_mutation_index(self)
        self.Assm_model.Assm_view.print_mutation_long_index(self)
        
    def addMutations(self):
        for gene_str in self.Graph:
            if gene_str not in self.mutations:
                continue # No predefined mutation in this gene
            for mutations_str in self.mutations[gene_str]:
                mu_type, coor, original, mutated = self.splitMutationInfo(gene_str, mutations_str)
                self.addNodeAndIdex(gene_str, mu_type, coor, original, mutated)
                    
    def splitMutationInfo(self, gene_str, mutations_str):
        mu_type = self.mutations[gene_str][mutations_str]
        mutation_tuple = mutations_str.split(':')
        original = mutation_tuple[0]
        mutated = mutation_tuple[1]
        #chr_str = mutation_tuple[2]
        try:
            coor = int(mutation_tuple[3])
        except:
            print "Invalid mutation information:", mutations_str
        return mu_type, coor, original, mutated
        
    def addNodeAndIdex(self, gene_str, mu_type, coor, original, mutated):
        ## Adding new nodes for mutations ##
        new_node_idx = len(self.Graph[gene_str].G)
        self.Graph[gene_str].addNode(new_node_idx, mutated, self.default_cov, mu_type)
        if mu_type == 'INS':
            self.Graph[gene_str].addEdge(coor-1, new_node_idx, self.default_ref_cov)
            self.Graph[gene_str].addEdge(new_node_idx, coor, self.default_ref_cov)
        elif mu_type == 'DEL':
            self.Graph[gene_str].addEdge(coor-1, new_node_idx, self.default_ref_cov)
            self.Graph[gene_str].addEdge(new_node_idx, coor+len(original), self.default_ref_cov)
        elif mu_type == 'SUB' or mu_type == 'SNP':
            self.Graph[gene_str].addEdge(coor-1, new_node_idx, self.default_ref_cov)
            self.Graph[gene_str].addEdge(new_node_idx, coor+len(mutated), self.default_ref_cov)
        else:
            print "Can't parse mutation format:", gene_str, mu_type, coor, original, mutated
            return False        
        ## Adding USR indexes ##
        if mu_type == 'DEL': # For DEL, we pass the original string since the new string is empty
            self.buildUSRindex(gene_str, mu_type, original, coor, new_node_idx)
        elif len(mutated) >= 5: # USR index for long mutations
            self.buildLongUSRindex(gene_str, mu_type, mutated, coor, new_node_idx)
        else: # USR index for SNPs and short mutations
            self.buildUSRindex(gene_str, mu_type, mutated, coor, new_node_idx)
        return True

        
    ## Building indexes for short mutations (< 5bp) ##
    def buildUSRindex(self, gene_str, mu_type, mutated_str, coor, new_node_idx):
        mu_len = len(mutated_str)
        if mu_type == 'SNP' or mu_type == 'SUB':
            next_coor = coor + mu_len #+ 1
            for rel_start_idx in xrange(6-mu_len+1):
                key_seq, idx_list = self.getSixmerIdxesSUB(gene_str, coor, next_coor, rel_start_idx, new_node_idx, mutated_str)                
                self.addSixmerToIndex(gene_str, key_seq, idx_list)
        elif mu_type == 'INS':
            next_coor = coor #+ 1
            for rel_start_idx in xrange(6-mu_len+1):
                key_seq, idx_list = self.getSixmerIdxesINS(gene_str, coor, next_coor, rel_start_idx, new_node_idx, mutated_str)                
                self.addSixmerToIndex(gene_str, key_seq, idx_list)
        elif mu_type == 'DEL':
            next_coor = coor + mu_len
            for rel_start_idx in xrange(1,6):
                key_seq, idx_list = self.getSixmerIdxesDEL(gene_str, coor, next_coor, rel_start_idx, new_node_idx, mutated_str)                
                self.addSixmerToIndex(gene_str, key_seq, idx_list)
                
    def getSixmerIdxesSUB(self, gene_str, first_end, next_start, rel_start_idx, new_node_idx, mutated_str):
        target_seq = self.Assm_model.target_sequences[gene_str]
        key_seq = ''
        idx_list = []
        for i in xrange(first_end - rel_start_idx, first_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        idx_list.append(new_node_idx)
        key_seq += mutated_str
        next_end = next_start + 6 - len(key_seq)
        for i in xrange(next_start, next_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        if len(key_seq) != 6: # or len(idx_list) != 6:
            print "Debug: getSixmerIdxesSUB()!!!!", 'key_seq:', key_seq
        return key_seq, idx_list

    def getSixmerIdxesDEL(self, gene_str, first_end, next_start, rel_start_idx, new_node_idx, mutated_str):
        target_seq = self.Assm_model.target_sequences[gene_str]
        key_seq = ''
        idx_list = []
        for i in xrange(first_end - rel_start_idx, first_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        next_end = next_start + 6 - rel_start_idx
        for i in xrange(next_start, next_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        idx_list.append(new_node_idx) # only for deletions, new node goes to the end
        if len(key_seq) != 6  or len(idx_list) != 7:
            print "Debug: getSixmerIdxesDEL()!!!!", 'key_seq:', key_seq
        return key_seq, idx_list

    def getSixmerIdxesINS(self, gene_str, first_end, next_start, rel_start_idx, new_node_idx, mutated_str):
        target_seq = self.Assm_model.target_sequences[gene_str]
        key_seq = ''
        idx_list = []
        for i in xrange(first_end - rel_start_idx, first_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        idx_list.append(new_node_idx)
        key_seq += mutated_str
        next_end = next_start + 6 - len(key_seq)
        for i in xrange(next_start, next_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        if len(key_seq) != 6: # or len(idx_list) != 6:
            print "Debug: getSixmerIdxesINS()!!!!", 'key_seq:', key_seq
        return key_seq, idx_list

    def addSixmerToIndex(self, gene_str, seq, idx_list):
        if seq in self.sixmer_dic_mutation:
            if gene_str in self.sixmer_dic_mutation[seq]:
                self.sixmer_dic_mutation[seq][gene_str].append(idx_list)
            else:
                self.sixmer_dic_mutation[seq][gene_str] = []
                self.sixmer_dic_mutation[seq][gene_str].append(idx_list)
        else:
            self.sixmer_dic_mutation[seq] = {}
            self.sixmer_dic_mutation[seq][gene_str] = []
            self.sixmer_dic_mutation[seq][gene_str].append(idx_list)
            
    def Coverage(self, molecule_inst, gene_str):
        new_mutation_search_list = []
        sixmers_with_mu = molecule_inst.mutation_search_list
        for sixmer in sixmers_with_mu:
            if sixmer in self.sixmer_dic_mutation and \
               gene_str in self.sixmer_dic_mutation[sixmer]:
                    normalized_vote = 1.0 / len(self.sixmer_dic_mutation[sixmer][gene_str])
                    self.iterLocation(molecule_inst, gene_str, sixmer, normalized_vote)
                    molecule_inst.mutation_barcode_cnt += 1
                    molecule_inst.mutation_barcode_list.append(sixmer)
            else: # When there is no perfect match, append to mutation search list
                new_mutation_search_list.append(sixmer)
        molecule_inst.mutation_search_list = new_mutation_search_list
        
    def iterLocation(self, molecule_inst, gene_str, sixmer, normalized_vote):
        for pos_ind in xrange(len(self.sixmer_dic_mutation[sixmer][gene_str])):
            six_positions = self.sixmer_dic_mutation[sixmer][gene_str][pos_ind]
            for idx in six_positions: # Increase Coverages
                if idx < 0:
                    continue
                self.Assm_model.g_Graph[gene_str].G[idx].cov += normalized_vote #Update Gloabal Graph
                molecule_inst.Graph.G[idx].cov += normalized_vote #Increase molecule coverage
        
    def iterMaskedLocation(self, molecule_inst, gene_str, sixmer, normalized_vote, dont_increase):
        for pos_ind in xrange(len(self.sixmer_dic_mutation[sixmer][gene_str])):
            six_positions = self.sixmer_dic_mutation[sixmer][gene_str][pos_ind]
            for i in xrange(len(six_positions)): # Increase Coverages
                idx = six_positions[i]
                if i in dont_increase or idx < 0:
#                    if len(six_positions) == 7:
#                        print i, idx, dont_increase, six_positions
                    continue
                self.Assm_model.g_Graph[gene_str].G[idx].cov += normalized_vote #Update Gloabal Graph
                molecule_inst.Graph.G[idx].cov += normalized_vote #Increase molecule coverage
            
    def CoverageLong(self, molecule_inst, gene_str):
        new_mutation_search_list = []
        sixmers_with_mu = molecule_inst.mutation_search_list
        for sixmer in sixmers_with_mu:
            if sixmer in self.sixmer_dic_l_mutation and \
               gene_str in self.sixmer_dic_l_mutation[sixmer]:
                    normalized_vote = 1.0 / len(self.sixmer_dic_l_mutation[sixmer][gene_str])
                    self.iterLocationLong(molecule_inst, gene_str, sixmer, normalized_vote)
            else: # When there is no perfect match, append to mutation search list
                new_mutation_search_list.append(sixmer)
        molecule_inst.mutation_search_list = new_mutation_search_list
                
    def iterLocationLong(self, molecule_inst, gene_str, sixmer, normalized_vote):
        for pos_ind in xrange(len(self.sixmer_dic_l_mutation[sixmer][gene_str])):
            six_positions = self.sixmer_dic_l_mutation[sixmer][gene_str][pos_ind]
            for idx in six_positions: # Increase Coverages
                self.Assm_model.g_Graph[gene_str].G[idx].cov += normalized_vote #Update Gloabal Graph
                molecule_inst.Graph.G[idx].cov += normalized_vote #Increase molecule coverage 
                
    ## Building indexes for long mutations (>= 5bp) ##
    def buildLongUSRindex(self, gene_str, mu_type, mutated_str, coor, new_node_idx):
        mu_len = len(mutated_str)
        if mu_type == 'SNP' or mu_type == 'DEL':
            print "SNP and DEL must be handled in buildUSRindex(). Not in buildLongUSRindex()"
        elif mu_type == 'SUB':
            next_coor = coor + mu_len #+ 1
            for rel_start_idx in xrange(1,6):
                self.getSixmerIdxesLongMut(gene_str, coor, next_coor, rel_start_idx, new_node_idx, mutated_str)
        elif mu_type == 'INS':
            next_coor = coor #+ 1
            for rel_start_idx in xrange(1,6):
                self.getSixmerIdxesLongMut(gene_str, coor, next_coor, rel_start_idx, new_node_idx, mutated_str)
        ## Add long mutation (itself) seq indexes ##
        self.getSixmerIdxesLongMutItself(gene_str, coor, new_node_idx, mutated_str)

    def getSixmerIdxesLongMut(self, gene_str, first_end, next_start, rel_start_idx, new_node_idx, mutated_str):
        target_seq = self.Assm_model.target_sequences[gene_str]
        ## Add 5-prime indexes ##
        key_seq = ''
        idx_list = []
        for i in xrange(first_end - rel_start_idx, first_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        idx_list.append(new_node_idx)
        key_seq += mutated_str[0:6-rel_start_idx]
        if len(key_seq) != 6:
            print "Debug: getSixmerIdxesLongMut()!!!!", 'key_seq:', key_seq
        self.addLongMutationIndex(gene_str, key_seq, idx_list)
        ## Add 3-prime indexes ##
        key_seq = ''
        idx_list = []
        idx_list.append(new_node_idx)
        key_seq += mutated_str[-rel_start_idx:]
        next_end = next_start + 6 - rel_start_idx
        for i in xrange(next_start, next_end):
            key_seq += target_seq[i]
            idx_list.append(i)
        if len(key_seq) != 6:
            print "Debug: getSixmerIdxesLongMut()!!!!", 'key_seq:', key_seq
        self.addLongMutationIndex(gene_str, key_seq, idx_list)

    def getSixmerIdxesLongMutItself(self, gene_str, coor, new_node_idx, mutated_str):
        for i in xrange(0, len(mutated_str)-5):
            key_seq = ''
            idx_list = []
            key_seq = mutated_str[i:i+6]
            idx_list.append(new_node_idx)
            self.addLongMutationIndex(gene_str, key_seq, idx_list)
    
    def addLongMutationIndex(self, gene_str, seq, idx_list):
        if seq in self.sixmer_dic_l_mutation:
            if gene_str in self.sixmer_dic_l_mutation[seq]:
                self.sixmer_dic_l_mutation[seq][gene_str].append(idx_list)
            else:
                self.sixmer_dic_l_mutation[seq][gene_str] = []
                self.sixmer_dic_l_mutation[seq][gene_str].append(idx_list)
        else:
            self.sixmer_dic_l_mutation[seq] = {}
            self.sixmer_dic_l_mutation[seq][gene_str] = []
            self.sixmer_dic_l_mutation[seq][gene_str].append(idx_list)
            
            
    ## Parsing mutations from predefined mutation file ##
    def parseMutationFile(self): # Testing/Debugging mutations replace it to parseMutationFile()
        with open(self.mutation_file_name, 'r') as inFile:
            for line in inFile:
                line = line.strip()
                if line.startswith('#') or line == '':
                    continue
                data = line.split(':')
                tmp_gene = data[0]
                original_str = data[1]
                mutated_str = data[2]
                chr_str = data[3]
                start_pos = data[4] # 0-based
                mutation_type = data[5]
                if tmp_gene not in self.mutations:
                    self.mutations[tmp_gene] = {}
                tmp_key = original_str + ':' + mutated_str + ':' + chr_str + ':' + start_pos
                self.mutations[tmp_gene][tmp_key] = mutation_type

    ## Print for testing
    def printGraphTest(self, gene):
        Single_Assm_model = self.Graph[gene].G
        out_file_name = os.path.join(self.Assm_model.Assm_view.path.out_dir, 'out_G_test.txt')
        with open(out_file_name,'w') as oFile: # print graph structure
            oFile.write('>\n')
            for idx, node in Single_Assm_model.iteritems():
                #print node.idx, node.base, node.from_edges, node.to_edges, node.cov
                oFile.write(str(node.idx) + '\t' +
                    str(node.base) + '\t' +
                    'from:' +
                    str(node.from_edges) + '\t' +
                    'to:' +
                    str(node.to_edges) + '\t' +
                    str(node.type) + '\t' +
                    str(node.cov) + '\n' )
                
if __name__=='__main__':
    print "no main function defined for Mutation_modules.py"



# Mutation file example
#BRAF	NM_004333	1123	c.1784T>C	p.F595S	-	7:140453151-140453151	CHP2_BRAF_2	TCAAACTGATGGGACCCACTCCATCGAGATTTC
#BRAF	NM_004333	1126	c.1789_1790CT>TC	p.L597S	-	7:140453145-140453146	CHP2_BRAF_2	TCAAACTGATGGGACCCACTCCATCGA
#BRAF	NM_004333	1127	c.1797_1799AGT>GAG	p.V600R	-	7:140453136-140453138	CHP2_BRAF_2	TCAAACTGATGGGACCCACTC
#BRAF	NM_004333	1128	c.1797_1797A>TACTACG	p.T599_V600insTT	-	7:140453138-140453138	CHP2_BRAF_2	TCAAACTGATGGGAC
#BRAF	NM_004333	1133	c.1799_1801delTGA	p.V600_K601>E	-	7:140453134-140453136	CHP2_BRAF_2	TCAAACTGATGGGACCCACTC
#BRAF	NM_004333	1135	c.1813_1814AG>TT	p.S605F	-	7:140453121-140453122	CHP2_BRAF_2	TCAAACTGATGGGACCCACTCCATCGA
#BRAF	NM_004333	144982	c.1797_1798insACA	p.T599_V600insT	-	7:140453137-140453138	CHP2_BRAF_2	TCAAACTGATGGGAC

# MAF example
#KRAS	3845	broad.mit.edu	37	12	25398284	25398285	+	
#Missense_Mutation	DNP	CC	CC	AA	rs121913529		
#TCGA-78-7167-01A-11D-2063-08	TCGA-78-7167-11A-01D-2063-08	CC	
#CC	-	-	CC	CC	Unknown	Untested	Somatic	
