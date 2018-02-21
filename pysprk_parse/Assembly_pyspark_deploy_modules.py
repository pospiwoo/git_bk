# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:30:24 2017

@author: Sunghee Woo
"""
import os#, sys#, time, operator, logging, copy #
import glob#, random, tempfile
#from operator import add
from collections import defaultdict
import shutil
import Assembly_modules as Assm
import Parameters as param
import Path_modules as Path
#from pyspark.sql import Row
from pyspark.sql import SQLContext
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
path = Path.Path()


class AssemblyPysparkProcessing:
    def __init__(self, fastq, target_fa):
        self.fastq = fastq
        self.target_fa = target_fa
        self.target_cnt_dic = self.parseTargetList()
        self.target_len_dic = self.parseTargetList()
        Assm_view = Assm.AssemblyView(path)
        Assm_model = Assm.AssemblyModel(fastq, target_fa, Assm_view)
        self.Assm_controller = Assm.AssemblyController(Assm_model, Assm_view)
        self.pos_str = ''
        self.newline_str = '\n'
        self.newline_str = '*'
        output_path_tmp = os.path.basename(self.fastq).split('.')[0]
        #self.output_path = '/home/swoo/AWS/output_pyspark_testing/'+output_path_tmp+'_'+str(os.getpid())
        #self.output_path = 's3://shortstackdeploy/'+output_path_tmp+'_'+str(os.getpid())
        self.output_path = 's3://shortstackdeploy/Output'
    
    def parseHeader(self):
        with open(self.fastq,"r") as inFile:
            first_header_line = inFile.readline()
            self.Assm_controller.Assm_model.ParseModel.parseColorEncodingTable(self.Assm_controller.Assm_model.color_encoding_file_name)            
            data = first_header_line.strip().strip(',')
            data = data.split(',')
            self.Assm_controller.Assm_model.ParseModel.idxToPoolList = [-1] * len(data)
            for i in xrange(len(data)):
                if 'P' in data[i]:
                    pool = data[i].split('P')[-1].strip()
                    self.Assm_controller.Assm_model.ParseModel.idxToPoolList[i] = pool
    
    def processHeader(self):
        with open(self.fastq,"r") as inFile:
            first_header_line = inFile.readline()
            self.Assm_controller.Assm_model.ParseModel.parseColorEncodingTable(self.Assm_controller.Assm_model.color_encoding_file_name)            
            data = first_header_line.strip().strip(',')
            data = data.split(',')
            self.Assm_controller.Assm_model.ParseModel.idxToPoolList = [-1] * len(data)
            for i in xrange(len(data)):
                if 'P' in data[i]:
                    pool = data[i].split('P')[-1].strip()
                    self.Assm_controller.Assm_model.ParseModel.idxToPoolList[i] = pool

    def parseTargetList(self):
        tmp_dic = {}
        with open(self.target_fa,'r') as targetFASTA:
            for line in targetFASTA:
                if line.startswith('>'):
                    curr_gene = line.strip()[1:]
                    tmp_dic[curr_gene] = ''
                else:
                    #tmp_dic[curr_gene] = [0.0]*len(line.strip())
                    tmp_dic[curr_gene] = len(line.strip())
        return tmp_dic
        
    def putIntoMolecule(self, line, mole_model):        
        data = line.strip().split(',')
        feat = data[0].strip()
        x = data[1].strip()
        y = data[2].strip()
        self.pos_str = x + '.' + y + '_' + feat
#        print '\n\n\n\n\n\n\n', len(data), len(self.idxToPoolList),'\n\n\n\n\n\n\n'
        for i in xrange(3, len(data)):
            pool = self.Assm_controller.Assm_model.ParseModel.idxToPoolList[i]
            if data[i] == '0': # No information
                continue
            elif len(data[i]) == 4:
                str_to_decode = '00' + data[i]
            elif len(data[i]) == 6:
                str_to_decode = data[i]
            else: # ignore half colors and less than 2
                continue
            self.Assm_controller.Assm_model.ParseModel.decode(pool, str_to_decode, mole_model)
    
        
    def putIntoMoleculeTSV(self, line, mole_model):        
        data = line.strip().split('\t')
        feat = data[0].strip()
        x = data[1].strip()
        y = data[2].strip()
        self.pos_str = x + '.' + y + '_' + feat
        
        for i in xrange(3, len(data)):
            if data[i].count('N') > 2:
                continue
            elif data[i].find('NN') > -1:
                mole_model.sixmer_NN_list.append(data[i])
                continue
            else:
                mole_model.sixmer_list.append(data[i])
        self.cnt_filtered_out_sixmers = 0
        self.cnt_used_sixmers = 0
        
        
    def genFASTAstr(self, Graph_inst, header_str, gene_str):
        #return_list = []
        return_str = ''
        if header_str:
            return_str += '>'+gene_str+';'+header_str+self.newline_str
        else:
            return_str += '>'+gene_str+self.newline_str
        return_str += ''.join(Graph_inst.max_seq)+self.newline_str
        return_str += ''.join(Graph_inst.max_seq_types)+self.newline_str
        if Graph_inst.max_seq_score != []:
            return_str += ''.join(Graph_inst.max_seq_score)+self.newline_str
            
#        for i in xrange(0,len(self.target_cnt_dic)):
    
        return return_str
    
    def genCOVstr(self, Graph_inst):
        return_str = ''
        for idx, node in Graph_inst.G.iteritems():
            if node.type == 'ref':
                return_str += str(node.cov) + self.newline_str
        return return_str.strip(self.newline_str)
#        return_list = []
#        if Graph_inst.max_seq_score != []:
#            return_list = Graph_inst.max_seq_cov
#        for idx, node in Graph_inst.G.iteritems():
#            if node.type == 'ref':
#                return_list.append(str(node.cov))
#        return return_list
    
    def genCOVlist(self, Graph_inst):
        return_list = []
        for idx, node in Graph_inst.G.iteritems():
            if node.type == 'ref':
                return_list.append(str(node.cov))
        return return_list
            
    def parseAndProcessCSVPyspark(self, line):
        if line.startswith('Features'):
            return '' # if header line, return
        molecule_model = Assm.MoleculeModel(self.Assm_controller.Assm_model)
        
        if self.fastq.endswith('.csv'):
            self.putIntoMolecule(line, molecule_model)
        elif self.fastq.endswith('.tsv'):
            self.putIntoMoleculeTSV(line, molecule_model)
        else: # Unknown file format
            return
            
        if self.Assm_controller.Assm_model.ParseModel.checkInsufficientSixmers(molecule_model):
            return ''
        
        origin_gene, MTM_ambig = self.Assm_controller.determineTargetGene(molecule_model)
        if MTM_ambig == True:
            self.Assm_controller.MTM_removed_cnt += 1
            return ''
        elif origin_gene.startswith('XXX') and not param.show_XXX_targets():
            self.Assm_controller.XXX_removed_cnt += 1
            return ''
        elif self.Assm_controller.Assm_model.ParseModel.checkInsufficientTargets(self.Assm_controller.Assm_model, molecule_model, origin_gene):
            return '' # Check if we have enough candidate targets
        
        molecule_model.initMutationGraph(origin_gene) # Charlie's Mutation-Graph-factory
        
        self.Assm_controller.estimateAllCov(molecule_model, origin_gene)
#        self.updateReadCounts(molecule_model, position) ### Update info ###
        
        if param.enable_blind_mu():
            self.Assm_controller.FindMutations(molecule_model, origin_gene)               
            molecule_model.Graph.GraphTrimming(self.Assm_controller.Assm_model.trimming_threshold)
            molecule_model.Graph.CallInsertions(self.Assm_controller.Assm_model.trimming_threshold)        
        
        if param.fast_path():
            molecule_model.Graph.grdQualityPath(0)
        else:
            molecule_model.Graph.optPath()            
        
        ### Update info ###
        position = ''
        self.updateInfo(molecule_model, position)
        
        FASTA_str = self.genFASTAstr(molecule_model.Graph, self.pos_str ,origin_gene)
        cov_list = self.genCOVstr(molecule_model.Graph)
#        cov_list = self.genCOVlist(molecule_model.Graph)
        VCF_list = ''
        VCF_list = ['']
        del molecule_model ### free molecule_model instance ###
#        cov_list = [int(1000*random.random()) for i in xrange(1)]
#        VCF_list = [int(1000*random.random()) for i in xrange(5)]
#        tmp_list2 = [int(1000*random.random()) for i in xrange(5)]
        return FASTA_str, origin_gene, cov_list, VCF_list
    
    def sum_cov(self, line1, line2):
        return_list = []
        try:
            float(line2[0])
            for i in xrange(0,len(line2)):
                return_list.append(line1[i] + line2[i])
        except :
            for i in xrange(0,len(line2)): 
    #            return_list.append(float(line1[i].encode('ascii','ignore')) + \
                return_list.append(line1[i] + \
                                   float(line2[i].encode('ascii','ignore')))
#        if line2[0].replace('.','',1).isdigit():
#        else:
        return return_list
    
    def updateInfo(self, molecule_model, position):
        molecule_model.pos_header = position
        self.num_molecules += 1
        self.num_reads += len(molecule_model.sixmer_list) + len(molecule_model.sixmer_NN_list)
        self.cnt_sixmer += len(molecule_model.sixmer_list)
        self.cnt_no_match_sixmer += len(molecule_model.mutation_search_list)
        self.cnt_NN_mer += len(molecule_model.sixmer_NN_list)
        
    def initAccumulators(self, sc_inst):
        self.num_molecules = sc_inst.accumulator(0)
        self.num_reads = sc_inst.accumulator(0)
        self.cnt_sixmer = sc_inst.accumulator(0)
        self.cnt_no_match_sixmer = sc_inst.accumulator(0)
        self.cnt_NN_mer = sc_inst.accumulator(0)
        self.cnt_filtered_out_sixmers = sc_inst.accumulator(0)
        self.cnt_used_sixmers = sc_inst.accumulator(0)
#        for gene in self.target_cnt_dic:
#            self.target_cnt_dic[gene] = sc.accumulator(0)
        
    def writeAccumulators(self, sc_inst):
        self.write_str = ''
        self.write_str += "\nTotal # molecules: " + str(self.num_molecules)
        self.write_str += "\nTotal # reads: " + str(self.num_reads)
        self.write_str += '\nTotal 3-spotter sixmers: ' + str(self.cnt_sixmer)
        self.write_str += "\nAvg 3-spotter sixmers per molecule:{0:.2f} ".format( \
                        float(str(self.cnt_sixmer))/float(str(self.num_molecules)))
        self.write_str += '\n3-spotter sixmers found perfect match: ' + \
                        str(int(str(self.cnt_sixmer)) - int(str(self.cnt_no_match_sixmer)))
        self.write_str += '\n3-spotter sixmers with no perfect match: ' + str(self.cnt_no_match_sixmer)
        self.write_str += '\nTotal 2-spotters:' + str(self.cnt_NN_mer)
        self.write_str += "\nAvg 2-spotter per molecule:{0:.2f} ".format( \
                        float(str(self.cnt_NN_mer))/float(str(self.num_molecules)))
#        self.write_str += '\nFiltered out 3-spotters: ' + str( \ 
#                        self.Assm_model.ParseModel.cnt_filtered_out_sixmers)
#        self.write_str += str(self.cnt_filtered_out_sixmers)
#        self.write_str += str(self.cnt_used_sixmers)
    
    
#    @execution_time
    def Process(self, sc):
        #self.parseHeader()
        self.initAccumulators(sc)
        
        text_file = sc.textFile(self.fastq)
        #header = text_file.first()
        #print '111111111111', header
        #self.processHeader()
        file_data = text_file.map(self.parseAndProcessCSVPyspark).filter(lambda x: x != '')
#                                 .map(lambda rtrn: Row(FASTA = rtrn[0], g1 = rtrn[1], g2=rtrn[2], g3=rtrn[3]))
        #with open("/home/swoo/AWS/test_file_out_3.txt","a") as oFile:
        #    oFile.write(file_data.collect())
        
        sqlContext = SQLContext(sc)
        file_dataframe = sqlContext.createDataFrame(file_data, ['FASTA','gene','coverage','VCF'])
#        file_dataframe.select('FASTA','gene','coverage').write.csv(self.output_path)
#        self.writeAccumulators(sc)

                               #Slow in Odin #
        file_dataframe.collect()
#        file_dataframe.select('FASTA').repartition(1).write.text(self.output_path)
        file_dataframe.select('FASTA').write.text(self.output_path)
        file_dataframe.select('gene','coverage').rdd.map(lambda row: [row[0], [x[0] for x in row[1]]])\
                            .aggregateByKey([0.0]*200, \
                                            self.sum_cov, self.sum_cov) \
#                            .repartition(1) \
                            .saveAsTextFile(os.path.join(self.output_path,'coverage'))
#                            .saveAsTextFile(os.path.join(self.output_path,'coverage'))
                               #Slow in Odin #
        
        
    def mergeFASTA(self):
        part_files = glob.glob(os.path.join(self.output_path,"*.txt"))
        with open(os.path.join(self.output_path,"out_FASTA.txt"),'wb') as oneFile:
            for p_file in part_files:
                with open(p_file,'rb') as partFile:
                    shutil.copyfileobj(partFile, oneFile)
                    
                    
    def mergeCOV(self):
        part_files = glob.glob(os.path.join(self.output_path,'coverage',"part*"))
        with open(os.path.join(self.output_path,"out_coverage.txt"),'wb') as oneFile:
            for p_file in part_files:
                with open(p_file,'rb') as partFile:
                    shutil.copyfileobj(partFile, oneFile)


    def mergeFASTAandSumCov(self):
        gene_cov_sum_dic = {}
        gene_cnt_dic = defaultdict(float)
        part_files = glob.glob(os.path.join(self.output_path,"*.csv"))
        with open(os.path.join(self.output_path,"out_FASTA.fa"),'wb') as outFASTAFile:
            for p_file in part_files:
                with open(p_file,'rb') as partFile:
                    for line in partFile:
                        data = line.replace('"','').strip().split(',')
                        gene_str = data[1]
                        cov_list = data[2].split(self.newline_str)
                        gene_cnt_dic[gene_str] += 1.0
                        if gene_str in gene_cov_sum_dic:
                            for i in xrange(0,len(cov_list)):
                                gene_cov_sum_dic[gene_str][i] += float(cov_list[i])
                        else:
                            gene_cov_sum_dic[gene_str] = []
                            for i in xrange(0,len(cov_list)):
                                gene_cov_sum_dic[gene_str].append(float(cov_list[i]))
                        outFASTAFile.write(data[0].replace(self.newline_str,'\n')+'\n')
                
        tmp_str = ''
        for gene in gene_cnt_dic:
            with open(os.path.join(self.output_path,"out_cov_avg_"+gene+".txt"),'wb') as outCovFile:
                tmp_write_str = gene+':'+str(int(gene_cnt_dic[gene]))+'\n'
                tmp_str += tmp_write_str
                newList = []
                for i in xrange(0,len(gene_cov_sum_dic[gene])):
                    tmp_write_str += str(gene_cov_sum_dic[gene][i]/gene_cnt_dic[gene]) + '\n'
                    newList.append(gene_cov_sum_dic[gene][i]/gene_cnt_dic[gene])
                outCovFile.write(tmp_write_str)
                
                plt.plot(newList)
                plt.ylabel('average coverage')
                plt.savefig(os.path.join(self.output_path,"out_cov_avg_"+gene+".pdf"))
                plt.close()

        self.write_str = tmp_str + '\n' + self.write_str


    def writeInfo(self):
        with open(os.path.join(self.output_path,"log.txt"),'wb') as outputMessageFile:
            outputMessageFile.write(self.write_str)




