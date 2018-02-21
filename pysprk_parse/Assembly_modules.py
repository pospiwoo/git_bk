# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 09:02:03 2016

@author: Sunghee Woo
"""
import os, time, operator, logging, copy #, sys
from datetime import timedelta
from collections import defaultdict
from scipy import stats
import numpy as np
import Alignment_modules as Align
import Graph_modules as Graph
import Mutation_modules as Mutation
import Parse_modules as Parse
import Parameters as param
import Path_modules as Path
path = Path.Path()

############### Model classes #################
# Molecule structure (Local)
class MoleculeModel:
    def __init__(self, Assm_model_inst):
        self.sixmer_dic_mole = {}
        self.sixmer_list = []
        self.sixmer_N_list = []
        self.sixmer_NN_list = []
        self.mutation_search_list = []
        self.Assm_model = Assm_model_inst
        self.Graph = {}
        self.gene_str = ''
        self.default_ref_cov = param.default_ref_cov()
        self.pos_header = ''
        self.cnt_reads = 0
        self.gene_uniq_sixmer = 0
        self.mutation_barcode_cnt = 0
        self.mutation_barcode_list = []

    def initGraph(self, gene_str):
        self.gene_str = gene_str
        seq_to_search = self.Assm_model.target_sequences[gene_str]
        for i in xrange(len(seq_to_search)):
            self.Graph.addNode(i, seq_to_search[i], self.default_ref_cov, 'ref')
        for i in xrange(1, len(seq_to_search)-1):
            self.Graph.addEdge(i-1, i, self.default_ref_cov)
            self.Graph.addEdge(i, i+1, self.default_ref_cov)

    def initMutationGraph(self, gene_str):
        self.gene_str = gene_str        
        self.Graph = copy.deepcopy(self.Assm_model.MutationModel.Graph[gene_str])

    def addSNP(self, position, snp_char, cov):
        from_node = position - 1
        to_node = position + 1
        original_seq_len = len(self.Assm_model.target_sequences[self.gene_str])
        self.Graph.lookupMutationNode(original_seq_len, from_node, to_node, snp_char, cov, 'SNP')

# Gloabal Assembler structure (Global)
class AssemblyModel:
    def __init__(self, fastq, target_fa, Assm_view_inst):
        self.num_molecules = 0
        self.num_reads = 0
        self.default_ref_cov = param.default_ref_cov()
        self.trimming_threshold = param.trimming_threshold()
        self.Assm_view = Assm_view_inst
        self.mutation_file_name = param.mutation_file_name()
        self.color_encoding_file_name = param.color_encoding_file_name()
        self.fastq_name = fastq
        self.ParseModel = Parse.ParseModel()
        self.is_FASTQ = self.ParseModel.checkInputFileType(fastq)
        self.target_fa_name = target_fa
        self.sixmer_dic_target = {}
        self.target_sequences = {}
        self.g_Graph = {}
        self.found_target_genes = defaultdict(int)
        self.Align_module = Align.Alignment()
        self.buildUSRindex()
        Graph_view = Graph.GraphView(self.Assm_view.path)
        self.MutationModel = Mutation.MutationModel(self, Graph_view)
        self.buildGlobalMutationGraph()  # Initiate global Mutation Graph
        self.ideal_barcode_cov = {}
        
#        logging.info("--- Start ---")
#        logging.info("FASTQ file name: " + fastq)
#        logging.info("target FASTA file name: " + target_fa)
        
    def buildUSRindex(self):
        #>BRAF_1
        #TTGTAGACTGTTCCAAATGATCCAGATCCAATTCTTTGTCCCACTGTAATCTGCCCATCAGGAATCT...
        with open(self.target_fa_name, mode='r') as fa:
            for line in fa:
                if line.find('>')>-1:
                    gene = line.strip()
                    gene_str = gene[1:]
                else:
                    seq = line.strip()
                    self.target_sequences[gene_str] = seq # Store all target sequences
                    for i in xrange(len(seq)-5):
                        key_seq = seq[i:i+6]
                        if key_seq in self.sixmer_dic_target:
                            if gene_str in self.sixmer_dic_target[key_seq]:
                                self.sixmer_dic_target[key_seq][gene_str].append(i)
                            else:
                                self.sixmer_dic_target[key_seq][gene_str] = []
                                self.sixmer_dic_target[key_seq][gene_str].append(i)
                        else:
                            self.sixmer_dic_target[key_seq] = {}
                            self.sixmer_dic_target[key_seq][gene_str] = []
                            self.sixmer_dic_target[key_seq][gene_str].append(i)
            self.Assm_view.print_target_index(self)

    def initGlobalGraph(self, target_seq):
        Graph_view = Graph.GraphView(self.Assm_view.path.out_dir)
        G_graph_per_gene = Graph.SpliceGraph(Graph_view)
        for i, seq in enumerate(target_seq):
            G_graph_per_gene.addNode(i, seq, self.default_ref_cov, 'ref')
        for i in xrange(1, len(target_seq)-1):
            G_graph_per_gene.addEdge(i-1, i, self.default_ref_cov)
            G_graph_per_gene.addEdge(i, i+1, self.default_ref_cov)
        return G_graph_per_gene

    def addGlobalSNP(self, gene_str, position, snp_char, cov):
        from_node = position - 1
        to_node = position + 1
        original_seq_len = len(self.target_sequences[gene_str])
        self.g_Graph[gene_str].lookupMutationNode(original_seq_len, from_node, to_node, snp_char, cov, 'SNP')        
        
    def buildGlobalMutationGraph(self):
        for gene_str in self.target_sequences:
            self.g_Graph[gene_str] = copy.deepcopy(self.MutationModel.Graph[gene_str])

############### Controller class #################
class AssemblyController:
    def __init__(self, Assm_model_inst, Assm_view_inst):
        self.Assm_model = Assm_model_inst
        self.Assm_view = Assm_view_inst
        self.NN_subst = ['AA', 'AT', 'AC', 'AG',\
                         'TA', 'TT', 'TC', 'TG',\
                         'CA', 'CT', 'CC', 'CG',\
                         'GA', 'GT', 'GC', 'GG']
        self.cnt_sixmer = 0
        self.cnt_no_match_sixmer = 0
        self.cnt_NN_mer = 0
        self.sum_zscores = []
        self.cnt_zscores = 0
        self.delta_zscore = 0.0
        self.MTM_removed_cnt = 0
        self.XXX_removed_cnt = 0
        self.MTM_score = {}
        for gene_str in self.Assm_model.target_sequences:
            self.MTM_score[gene_str] = 0.0
            self.sum_zscores.append(0.0)
        if param.debug_output_on():
            self.debugFile = open(param.get_debug_out_file_name(),'w')

    def execution_time(func):
        def wrapper(*args, **kwargs):
            begin = time.time()
            func(*args, **kwargs)
            finish = time.time()
            time_in_sec = finish - begin
            print 'Excecution time:', str(timedelta(seconds=time_in_sec))
        return wrapper
    
    @execution_time
    def Process(self):
        if param.ideal_cov():
            self.processIdeal()
        
        with open(self.Assm_model.fastq_name,'r') as fq, \
            open(os.path.join(self.Assm_view.path.out_dir, param.assembled_read_out()),'w') as oFile:
            
            ### Parse header (along with the first line) ###
            if self.Assm_model.is_FASTQ:
                flag_keep_reading, position, offset_track, headers = self.Assm_model.ParseModel.parseFQHeaders(fq)
            else: # If it's not FASTQ, we assume it is an internal format from image-processing
                flag_keep_reading = self.Assm_model.ParseModel.parseCSVHeaders(fq, self.Assm_model.color_encoding_file_name)
            
            while flag_keep_reading: # This block iterates per molecule block
                ### allocate molecule_model instance ###
                molecule_model = MoleculeModel(self.Assm_model)
                
                ### Parse each molecule at a time ###
                if self.Assm_model.is_FASTQ:
                    flag_keep_reading, position, offset_track = self.Assm_model.ParseModel.parseFQFile(fq, molecule_model, position, offset_track)
                else:
                    flag_keep_reading, position = self.Assm_model.ParseModel.parseCSVFile(fq, molecule_model)
                    if self.Assm_model.ParseModel.checkInsufficientSixmers(molecule_model):
                        continue # Check if we have enough sixmers
                
                ### Find perfect match and determine the originated gene target
                origin_gene, MTM_ambig = self.determineTargetGene(molecule_model)
                if MTM_ambig == True:
                    self.MTM_removed_cnt += 1
                    continue
                elif origin_gene.startswith('XXX') and not param.show_XXX_targets():
                    self.XXX_removed_cnt += 1
                    continue
                elif self.Assm_model.ParseModel.checkInsufficientTargets(self.Assm_model, molecule_model, origin_gene):
                    continue # Check if we have enough candidate targets
                
                ### Initiate Graph structure ###
                molecule_model.initMutationGraph(origin_gene) # Charlie's Mutation-Graph-factory
                
                ### Estimating coverage and Coloring the Graph ###
                self.estimateAllCov(molecule_model, origin_gene)
                
                ### Update info ###
                self.updateReadCounts(molecule_model, position)
                
                if param.enable_blind_mu():
                    ### Find mutations and add to mutation Graph ###
                    self.FindMutations(molecule_model, origin_gene)                    
                    ### Graph trimming ###
                    molecule_model.Graph.GraphTrimming(self.Assm_model.trimming_threshold)                    
                    ### Call insertions before computing the heaviest/longest path ###
                    molecule_model.Graph.CallInsertions(self.Assm_model.trimming_threshold)
                
                ### Find maximum scoring path ###
                if param.fast_path():
                    molecule_model.Graph.grdQualityPath(0)
                else:
                    molecule_model.Graph.optPath()
                
                molecule_model.Graph.GraphView.PrintFullReads(molecule_model.Graph, oFile, molecule_model.pos_header 
                                                              + ';3Spotters_cnt:' + str(len(molecule_model.sixmer_list))
                                                              + ';gene_uniq_sixmer:' + str(molecule_model.gene_uniq_sixmer) 
                                                              + ';mutation_search_cnt:' + str(len(molecule_model.mutation_search_list)) 
                                                              + ';2Spotters_cnt:' + str(len(molecule_model.sixmer_NN_list)) 
                                                              #+ ';mu_bar_cnt:' + str(molecule_model.mutation_barcode_cnt) 
                                                              + ';3Spotters:' + ','.join(molecule_model.sixmer_list)
                                                              + ';2Spotters:' + ','.join(molecule_model.sixmer_NN_list)
                                                              + ';Mutation3Spotters:' + ','.join(molecule_model.mutation_search_list)
                                                              ,origin_gene)
                #molecule_model.Graph.GraphView.PrintFullReadsDynamicP(molecule_model.Graph, oFile, molecule_model.pos_header, origin_gene)
                
                ### free molecule_model instance ###
                del molecule_model
            
            ###  Print all outputs ###
            self.DeltaZScore()
            self.Assm_view.print_outputs(self.Assm_model, self)
            print 'MTM zscore removed:', self.MTM_removed_cnt
            print 'Feature contaminants removed:', self.XXX_removed_cnt
            self.Assm_view.out_message_str += '\nMTM zscore removed: ' + str(self.MTM_removed_cnt)
            self.Assm_view.out_message_str += '\nFeature contaminants removed: ' + str(self.XXX_removed_cnt)
            with open(param.output_message(),'w') as outputMessageFile:
                outputMessageFile.write(self.Assm_view.out_message_str)
                
        if param.debug_output_on():
            self.debugFile.close()
            
    def processIdeal(self):
        ideal_file_name = param.ideal_file_name()
        with open(ideal_file_name,'r') as fq:
            flag_keep_reading = True
            while flag_keep_reading:
                molecule_model = MoleculeModel(self.Assm_model)      
                ### Parse each molecule at a time ###
                flag_keep_reading, origin_gene = self.Assm_model.ParseModel.parseBarcodeFile(fq, molecule_model)
                if origin_gene == '' or origin_gene.startswith('COSM'):
                    continue
                #self.Assm_model.found_target_genes[origin_gene] = 1
#                if not origin_gene.startswith('BRAF'):
#                    continue
                ### Initiate Graph structure ###
                molecule_model.initMutationGraph(origin_gene) # Charlie's Mutation-Graph-factory                
                ### Estimating coverage and Coloring the Graph ###
                self.estimateAllCov(molecule_model, origin_gene)                
                ### Find maximum scoring path ###
                molecule_model.Graph.greedyQualityPath(0)                
                self.Assm_model.ideal_barcode_cov[origin_gene] = molecule_model.Graph.max_seq_cov
                ideal_out_file_name = os.path.join(self.Assm_view.path.out_dir, origin_gene+'_ideal.txt')
                if param.print_ideal_cov():
                    idealFile = open(ideal_out_file_name,'w')
                    for j in xrange(0, len(molecule_model.Graph.max_seq_cov)):
                        idealFile.write(str(float(molecule_model.Graph.max_seq_cov[j])-self.Assm_model.default_ref_cov)+'\n')
                    idealFile.close()
            
    def estimateAllCov(self, molecule_model, origin_gene):
        self.Coverage(molecule_model, origin_gene) # Perfect match
        self.Assm_model.MutationModel.Coverage(molecule_model, origin_gene) # Predefined short mutations
        self.Assm_model.MutationModel.CoverageLong(molecule_model, origin_gene) # Predefined Long mutations
        self.MaskedCoverage(molecule_model, origin_gene) # Reads with 'N'
    
    def Coverage(self, molecule_model_inst, gene_key):
        for tmp_sixmer in molecule_model_inst.sixmer_list:
            if tmp_sixmer in self.Assm_model.sixmer_dic_target and \
               gene_key in self.Assm_model.sixmer_dic_target[tmp_sixmer]:
                    normalized_vote = 1.0 / len(self.Assm_model.sixmer_dic_target[tmp_sixmer][gene_key])
                    self.iterLocation(molecule_model_inst, gene_key, tmp_sixmer, normalized_vote, [])
            else: # When there is no perfect match, append to mutation search list
                molecule_model_inst.mutation_search_list.append(tmp_sixmer)
    
    def MaskedCoverage(self, molecule_model_inst, gene_key):
        for tmp_sixmer_with_N in molecule_model_inst.sixmer_NN_list:
            dont_increase = tmp_sixmer_with_N.find('NN')
            dont_increase = [dont_increase, dont_increase + 1]
            if tmp_sixmer_with_N.count('NN') != 1:
                print "Reads with NN expected!", tmp_sixmer_with_N
                continue
            if param.two_spotter_use() == 0:
                normalized_vote = 0.0
            elif param.two_spotter_use() == 1: # divide by 16
                normalized_vote = 1.0 / len(self.NN_subst)
            elif param.two_spotter_use() == 2: # divide by possible locations
                location_cnt = self.assessNormalization(molecule_model_inst, gene_key, tmp_sixmer_with_N)
                if location_cnt == 0:
                    continue
                normalized_vote = 1.0 / location_cnt
            for j in xrange(len(self.NN_subst)):
                tmp_sixmer = tmp_sixmer_with_N.replace('NN',self.NN_subst[j])
                # in sixmer_dic_target
                if tmp_sixmer in self.Assm_model.sixmer_dic_target and \
                gene_key in self.Assm_model.sixmer_dic_target[tmp_sixmer]:
                    if param.dont_use_2_spot_in_global():
                        self.iterLocation_2_spotters(molecule_model_inst, gene_key, tmp_sixmer, normalized_vote, dont_increase)
                    else:
                        self.iterLocation(molecule_model_inst, gene_key, tmp_sixmer, normalized_vote, dont_increase)
                # in sixmer_dic_mutation
                if tmp_sixmer in self.Assm_model.MutationModel.sixmer_dic_mutation and \
                gene_key in self.Assm_model.MutationModel.sixmer_dic_mutation[tmp_sixmer]:
                    self.Assm_model.MutationModel.iterMaskedLocation(molecule_model_inst, \
                    gene_key, tmp_sixmer, normalized_vote, dont_increase)
                    continue
                
    def assessNormalization(self, molecule_model_inst, gene_key, sixmer_with_N):
        accurence = 0
        for j in xrange(len(self.NN_subst)):
            tmp_sixmer = sixmer_with_N.replace('NN',self.NN_subst[j])
            # In sixmer_dic_target
            if tmp_sixmer in self.Assm_model.sixmer_dic_target and \
                gene_key in self.Assm_model.sixmer_dic_target[tmp_sixmer]:
                accurence += len(self.Assm_model.sixmer_dic_target[tmp_sixmer][gene_key])                
            # In sixmer_dic_mutation
            if tmp_sixmer in self.Assm_model.MutationModel.sixmer_dic_mutation and \
                gene_key in self.Assm_model.MutationModel.sixmer_dic_mutation[tmp_sixmer]:
                accurence += len(self.Assm_model.MutationModel.sixmer_dic_mutation[tmp_sixmer][gene_key])                
        return accurence
        
    def iterLocation(self, molecule_model_inst, gene_key, sixmer, normalized_vote, dont_increase):
        for pos_ind in xrange(len(self.Assm_model.sixmer_dic_target[sixmer][gene_key])):
            tmp_pos = self.Assm_model.sixmer_dic_target[sixmer][gene_key][pos_ind]
            # Increase Coverages
            self.increaseCov(molecule_model_inst, gene_key, tmp_pos, normalized_vote, dont_increase)
                
    def increaseCov(self, molecule_inst, gene_str, start_idx, normalized_vote, dont_increase):
        for i in xrange(0,6):
            if i in dont_increase:
                continue
            molecule_inst.Graph.G[start_idx+i].cov += normalized_vote #Increase molecule coverage
            self.Assm_model.g_Graph[gene_str].G[start_idx+i].cov += normalized_vote #Update Gloabal Graph
        
    def iterLocation_2_spotters(self, molecule_model_inst, gene_key, sixmer, normalized_vote, dont_increase):
        for pos_ind in xrange(len(self.Assm_model.sixmer_dic_target[sixmer][gene_key])):
            tmp_pos = self.Assm_model.sixmer_dic_target[sixmer][gene_key][pos_ind]
            # Increase Coverages
            self.increaseCov_2_spotters(molecule_model_inst, gene_key, tmp_pos, normalized_vote, dont_increase)
                
    def increaseCov_2_spotters(self, molecule_inst, gene_str, start_idx, normalized_vote, dont_increase):
        for i in xrange(0,6):
            if i in dont_increase:
                continue
            molecule_inst.Graph.G[start_idx+i].cov += normalized_vote #Increase molecule coverage
            
    def determineTargetGene(self, molecule_model_inst):
        MTM_score = copy.deepcopy(self.MTM_score)
        # map all sixmers to all targets
        for tmp_sixmer in molecule_model_inst.sixmer_list:
            if tmp_sixmer in self.Assm_model.sixmer_dic_target:
                for tmp_gene in self.Assm_model.sixmer_dic_target[tmp_sixmer]:
                    if param.is_genomic_dna():
                        MTM_score[tmp_gene] += len(self.Assm_model.sixmer_dic_target[tmp_sixmer][tmp_gene])
                    else:
                        MTM_score[tmp_gene] += 1.0
        # map all 2-spotters to all targets
        for tmp_sixmer_with_N in molecule_model_inst.sixmer_NN_list:
            for j in xrange(len(self.NN_subst)):
                tmp_sixmer = tmp_sixmer_with_N.replace('NN',self.NN_subst[j])
                if tmp_sixmer in self.Assm_model.sixmer_dic_target:
                    for tmp_gene in self.Assm_model.sixmer_dic_target[tmp_sixmer]:
                        MTM_score[tmp_gene] += 1.0/16
        MTM_score = sorted(MTM_score.items(), key=operator.itemgetter(1))
        if len(MTM_score) == 0: # Check if we got at leat one candidate gene
            return '', True
        genes = list(MTM_score)
        selected_gene = genes[-1][0]
        selected_gene_cnt = int(genes[-1][1])
        molecule_model_inst.gene_uniq_sixmer = selected_gene_cnt
        
        genes_a = []
        for j in xrange(0,len(genes)):
            genes_a.append(genes[j][1])
        genes_a = np.array(genes_a)
        zscore_list = stats.zscore(genes_a)
        
        if len(genes) > 1 and float(genes[-1][1]) - float(genes[-2][1]) < 0.1:
            return '', True
#        elif float(zscore_list[-1]) - float(zscore_list[-2]) < param.zscore_delta_threshold():
#            return '', True
        elif float(zscore_list[-1]) < float(param.zscore_threshold()):
            return '', True
        
        if not selected_gene.startswith('XXX'):
            self.cnt_zscores += 1.0
            for j in xrange(0,len(zscore_list)):
                if param.print_raw_z_score():
                    self.sum_zscores[j] += zscore_list[j]
                else:
                    self.sum_zscores[j] += genes[j][1]
        
        return selected_gene, False
        
    def FindMutations(self, molecule_inst, gene_str):
        sixmers_with_mu = molecule_inst.mutation_search_list
        target_seq = self.Assm_model.target_sequences[gene_str]
        for sixmer in sixmers_with_mu:
            snp_dic = self.Assm_model.Align_module.ProcessHamming(molecule_inst, target_seq, sixmer)
            if len(snp_dic) > 0: # SNPs found
                self.iterSNP(molecule_inst, gene_str, len(target_seq), snp_dic)
            else: # No SNPs found, find Insertions/Deletions
                #self.FindINS()
                #self.FindDel()
                continue

    def iterSNP(self, molecule_inst, gene_str, target_seq_len, snp_dic):
        normalized_vote = 1.0 / len(snp_dic)
        for rel_snp_idx in snp_dic:
            snp_char = snp_dic[rel_snp_idx][0]
            target_chnk = snp_dic[rel_snp_idx][1]
            locations = self.Assm_model.sixmer_dic_target[target_chnk][gene_str]
            for tmp_pos in locations:
                snp_idx = tmp_pos + rel_snp_idx
                if snp_idx <= 2 or snp_idx >= target_seq_len-3: #if snp_idx > 1 and snp_idx < len(target_seq)-1:
                    continue
                # Add SNP node and edges
                molecule_inst.addSNP(snp_idx, snp_char, normalized_vote)
                # Add Global SNP node and edges
                self.Assm_model.addGlobalSNP(gene_str, snp_idx, snp_char, normalized_vote)
                # Increase coverages
                self.increaseCov(molecule_inst, gene_str, tmp_pos, normalized_vote, [rel_snp_idx])
                
    def updateReadCounts(self, molecule_model, position):
        molecule_model.pos_header = position
        self.Assm_model.num_molecules += 1
        self.Assm_model.num_reads += molecule_model.cnt_reads
        self.cnt_sixmer += len(molecule_model.sixmer_list)
        self.cnt_NN_mer += len(molecule_model.sixmer_NN_list)
        self.cnt_no_match_sixmer += len(molecule_model.mutation_search_list)
#        print 'sixmers:', len(molecule_model.sixmer_list), 'NN_mers:', \
#                len(molecule_model.sixmer_NN_list), \
#                'no_match:', len(molecule_model.mutation_search_list), \
#                'match_found:', len(molecule_model.sixmer_list)-len(molecule_model.mutation_search_list)
#        print ''

    def DeltaZScore(self):
        for j in xrange(0,len(self.sum_zscores)):
            self.sum_zscores[j] = float(self.sum_zscores[j]) / self.cnt_zscores
        self.delta_zscore = self.sum_zscores[-1] - self.sum_zscores[-2]
        self.Assm_view.print_zscore(param.zscore_out(), self.sum_zscores)

############### View class #################
class AssemblyView:
    def __init__(self, path_inst):
        self.path = path_inst
        log_file_name = os.path.join(self.path.parent_dir, param.log_message())
#        logging.basicConfig(level=logging.DEBUG, filename=log_file_name, filemode='a+', \
#                        format='%(asctime)-15s %(levelname)-8s %(message)s')
        self.out_message_str = ''
            
    def print_matrix(self, matrix):
        for mat_line in matrix:
            print('\t'.join(map(str,mat_line)))
        return
        
    def print_target_index(self, model_inst):
        target_idx_file = open(param.USR_index_file() ,mode='w')
        for i in model_inst.sixmer_dic_target:
            target_idx_file.write(i+'\t'+str(len(model_inst.sixmer_dic_target[i]))+'\t')
            for j in model_inst.sixmer_dic_target[i]:
                target_idx_file.write(j+':')
                target_idx_file.write(','.join(map(str, model_inst.sixmer_dic_target[i][j])))
                target_idx_file.write('\t')
            target_idx_file.write('\n')
        target_idx_file.close()
        
    def print_mutation_index(self, model_inst):
        target_idx_file = open(param.USR_mutation_index_file() ,mode='w')
        for sixmer in model_inst.sixmer_dic_mutation:
            target_idx_file.write(sixmer+'\t'+str(len(model_inst.sixmer_dic_mutation[sixmer]))+'\t')
            for gene in model_inst.sixmer_dic_mutation[sixmer]:
                target_idx_file.write(gene+':')
                for loc in model_inst.sixmer_dic_mutation[sixmer][gene]:
                    target_idx_file.write('[')
                    target_idx_file.write(','.join(map(str, loc)))
                    target_idx_file.write(']')
                target_idx_file.write('\t')
            target_idx_file.write('\n')
        target_idx_file.close()
        
    def print_mutation_long_index(self, model_inst):
        target_idx_file = open(param.USR_long_mutation_index_file() ,mode='w')
        for sixmer in model_inst.sixmer_dic_l_mutation:
            target_idx_file.write(sixmer+'\t'+str(len(model_inst.sixmer_dic_l_mutation[sixmer]))+'\t')
            for gene in model_inst.sixmer_dic_l_mutation[sixmer]:
                target_idx_file.write(gene+':')
                print_str = ''
                for loc in model_inst.sixmer_dic_l_mutation[sixmer][gene]:
                    #target_idx_file.write(','.join(map(str, model_inst.sixmer_dic_l_mutation[sixmer][gene][loc])))
                    #target_idx_file.write(','.join(map(str, loc)))
                    for pos in xrange(len(loc)):
                        print_str += str(loc[pos]) + ','
                target_idx_file.write(print_str.strip(','))
                target_idx_file.write('\t')
            target_idx_file.write('\n')
        target_idx_file.close()
        
    def testParseByMolecule(self, molecule_model, tmp_position):
            print (tmp_position, molecule_model.sixmer_list, 
                len(molecule_model.sixmer_list), 
                len(molecule_model.sixmer_N_list), 
                len(molecule_model.sixmer_NN_list)
                )
            print tmp_position, molecule_model.sixmer_NN_list
        
    def print_time(self, start, end):
        time_in_sec = end - start
        print_str = 'Excecution time: ' + str(timedelta(seconds=time_in_sec))
#        logging.info(print_str)
#        logging.info("--- END ---")
        return
        
    def print_cov(self, cov_file_name, tmp_cov):
        cov_file_name = os.path.join(self.path.out_dir, cov_file_name)
        with open(cov_file_name, mode='w') as cov_file:
            for cov_val in tmp_cov:
                cov_file.write(str(cov_val)+"\n")
    
    def print_output_graphs(self, Assm_model, origin_gene):
        print origin_gene, ':', Assm_model.found_target_genes[origin_gene]
        self.out_message_str += '\n' + origin_gene + ': ' + str(Assm_model.found_target_genes[origin_gene])
        Single_Assm_model = Assm_model.g_Graph[origin_gene]
        Single_Assm_model.GraphView.PrintGraphMaxPath(Single_Assm_model, param.max_consensus_graph_out() + origin_gene + '.txt', origin_gene)
        #Single_Assm_model.GraphView.PrintGraph(Single_Assm_model, param.consensus_graph_out() + origin_gene + '.txt', origin_gene)
        if len(Assm_model.ideal_barcode_cov) > 0 and not origin_gene.startswith('XXX'):
            Single_Assm_model.GraphView.PrintGraphAvgWithIdeal(Single_Assm_model, param.avg_consensus_graph_out() + origin_gene + '.txt', Assm_model.found_target_genes[origin_gene], Assm_model.ideal_barcode_cov[origin_gene])
        else:
            Single_Assm_model.GraphView.PrintGraphAvg(Single_Assm_model, param.avg_consensus_graph_out() + origin_gene + '.txt', Assm_model.found_target_genes[origin_gene], origin_gene)
                
    def print_zscore(self, zscore_file_name, tmp_zscore):
        z_file_name = os.path.join(self.path.out_dir, zscore_file_name)
        z_delta_file_name = os.path.join(self.path.out_dir, 'delta_'+zscore_file_name)
        with open(z_file_name, mode='w') as z_file, \
            open(z_delta_file_name, mode='w') as z_d_file:
            for z in tmp_zscore:
                z_file.write(str(z)+"\n")
            for i in xrange(0,len(tmp_zscore)-1):
                z_d_file.write(str(tmp_zscore[i+1] - tmp_zscore[i])+"\n")
                
    def print_outputs(self, Assm_model, Assm_controller):
        for found_gene in Assm_model.found_target_genes:
            Assm_model.g_Graph[found_gene].greedyQualityPath(0)
            self.print_output_graphs(Assm_model, found_gene) # Output results to files
        print "Total # molecules", Assm_model.num_molecules
        print "Total # reads", Assm_model.num_reads
#        logging.info("Total # molecules: " + str(Assm_model.num_molecules))
#        logging.info("Total # reads: " + str(Assm_model.num_reads))
        print 'Total 3-spotter sixmers:', Assm_controller.cnt_sixmer
        print("Avg 3-spotter sixmers per molecule:{0:.2f}".format(float(Assm_controller.cnt_sixmer)/float(Assm_model.num_molecules)))
        print '3-spotter sixmers found perfect match:', Assm_controller.cnt_sixmer - Assm_controller.cnt_no_match_sixmer
        print("Avg perfect match 3-spotter sixmers per molecule:{0:.2f}".format(float(Assm_controller.cnt_sixmer - Assm_controller.cnt_no_match_sixmer)/float(Assm_model.num_molecules)))
        print '3-spotter sixmers with no perfect match:', Assm_controller.cnt_no_match_sixmer
        print 'Total 2-spotters:', Assm_controller.cnt_NN_mer
        print("Avg 2-spotter per molecule:{0:.2f}".format(float(Assm_controller.cnt_NN_mer)/float(Assm_model.num_molecules)))
        if param.filterout_orphans():
            print 'Filtered out 3-spotters:', Assm_controller.Assm_model.ParseModel.cnt_filtered_out_sixmers
            print 'Used 3-spotters:', Assm_controller.Assm_model.ParseModel.cnt_used_sixmers
        print 'Delta Z score:', Assm_controller.delta_zscore

        self.out_message_str += "\nTotal # molecules: " + str(Assm_model.num_molecules)
        self.out_message_str += "\nTotal # reads: " + str(Assm_model.num_reads)
        self.out_message_str += "\nTotal # molecules: " + str(Assm_model.num_molecules)
        self.out_message_str += "\nTotal # reads: " + str(Assm_model.num_reads)
        self.out_message_str += '\nTotal 3-spotter sixmers: ' + str(Assm_controller.cnt_sixmer)
        self.out_message_str += "\nAvg 3-spotter sixmers per molecule:{0:.2f} ".format(float(Assm_controller.cnt_sixmer)/float(Assm_model.num_molecules))
        self.out_message_str += '\n3-spotter sixmers found perfect match: ' + str(Assm_controller.cnt_sixmer - Assm_controller.cnt_no_match_sixmer)
        self.out_message_str += "\nAvg perfect match 3-spotter sixmers per molecule:{0:.2f} ".format(float(Assm_controller.cnt_sixmer - Assm_controller.cnt_no_match_sixmer)/float(Assm_model.num_molecules))
        self.out_message_str += '\n3-spotter sixmers with no perfect match: ' + str(Assm_controller.cnt_no_match_sixmer)
        self.out_message_str += '\nTotal 2-spotters:' + str(Assm_controller.cnt_NN_mer)
        self.out_message_str += "\nAvg 2-spotter per molecule:{0:.2f} ".format(float(Assm_controller.cnt_NN_mer)/float(Assm_model.num_molecules))
        if param.filterout_orphans():
            self.out_message_str += '\nFiltered out 3-spotters: ' + str(Assm_controller.Assm_model.ParseModel.cnt_filtered_out_sixmers)
            self.out_message_str += '\nUsed 3-spotters: ' + str(Assm_controller.Assm_model.ParseModel.cnt_used_sixmers)
        self.out_message_str += '\nDelta Z score: ' + str(Assm_controller.delta_zscore)

################ main #################
#if __name__=='__main__':
#    if len(sys.argv)<3:
#        print('Insufficient arguments')
#        print('Usage: python XXX.py <FASTQ> <Target_FASTA>')
#        exit()
#
#    fastq = sys.argv[1]
#    target_fa = sys.argv[2]
#
#    Assm_view = AssemblyView()
#    Assm_model = AssemblyModel(fastq, target_fa, Assm_view)
#    Assm_controller = AssemblyController(Assm_model, Assm_view)
#    
#    start = time.time()    
#    Assm_controller.Process() # Core module
#    end = time.time()
#    Assm_view.print_time(start, end)



