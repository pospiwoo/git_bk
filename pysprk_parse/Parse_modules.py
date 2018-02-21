# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 11:35:52 2016

@author: Sunghee Woo
"""
import os
import Path_modules as Path
import Parameters as param

class ParseModel:
    def __init__(self):
        self.encoding_dic = {}
        self.idxToPoolList = []      
        if param.filterout_orphans():
            if param.pool_spec_orphans():
                self.all_sixmer_file_name = 'copy_from_barcode_design_pool_spec.txt'
            else:
                self.all_sixmer_file_name = 'copy_from_barcode_design.txt'
            self.nonsense_sixmers = {}
            self.cnt_filtered_out_sixmers = 0
            self.all_sixmers = {}
            self.cnt_used_sixmers = 0
            self.path = Path.Path()
            self.all_sixmer_file_name = os.path.join(self.path.encoding_dir, self.all_sixmer_file_name)
            with open(self.all_sixmer_file_name,'r') as allsix_file:
                for line in allsix_file:
                    if param.pool_spec_orphans():
                        data = line.strip().split('\t')
                        tmp_pool = data[0]
                        tmp_sixmer = self.revcomp(data[1])
                        if not tmp_pool in self.all_sixmers:
                            self.all_sixmers[data[0]] = {}
                        self.all_sixmers[data[0]][tmp_sixmer] = 0
                    else:
                        self.all_sixmers[self.revcomp(line.strip())] = 0
            self.tmp_target_dic = {}

    ## From FASTQ ##
    def checkInputFileType(self, input_file_name):
        if input_file_name.endswith('.fastq') or input_file_name.endswith('.fq'):
            return True
        else:
            return False
        
    def parseFQHeaders(self, fq_file):
        headers = []
        for line in fq_file:
            line = line.strip()
            if line.startswith('#'):
                headers.append(line)
                continue
            else:
                header_str = line[1:]
                position_str = header_str.split('_')
                position_str = position_str[2]
                return True, position_str, 1, headers
        else: # No more lines to be read from file
            return False, '', -1, ''
        
    def parseFQFile(self, fq_file, molecule_model_inst, prev_position, offset):
        #@Seq_GridPos_568,2.1_Cycle=0_Pool=0
        #NNNNNN
        #+Seq_GridPos_568,2.1_Cycle=0_Pool=0
        #??????
        for line in fq_file:
            line = line.strip()
            if offset == 0:
                offset += 1
                header_str = line[1:]
                position_str = header_str.split('_')
                position_str = position_str[2]
    
                if position_str != prev_position:
                    return True, position_str, 1
                else:
                    continue
            elif offset == 1:
                offset += 1
                molecule_model_inst.cnt_reads += 1
                seq = line
                if seq.count('N') == 0:
                    molecule_model_inst.sixmer_list.append(seq)
                elif seq.count('N') == 1:
                    molecule_model_inst.sixmer_N_list.append(seq)
                elif seq.count('N') == 2:
                    molecule_model_inst.sixmer_NN_list.append(seq)
            elif offset == 2:
                offset += 1
                continue
            elif offset == 3:
                offset=0
                continue
            else:
                print 'Wrong FASTQ input', line
        else: # No more lines to be read from file
            return False, '', offset
            
    ## From image processed format ##
    def revcomp(self, seq):
        complement= {}
        complement['A']='T'
        complement['T']='A'
        complement['C']='G'
        complement['G']='C'
        complement['-']='-'
        complement['N']='N'
        revc=''
        if len(seq)==1:
            revc = complement[seq]
        else:
            for i in range(len(seq)):
                revc =  complement[seq[i]] + revc
        return(revc)
        
    def revcomp_dummy(self, seq):
        complement= {}
        complement['A']='T'
        complement['T']='A'
        complement['C']='G'
        complement['G']='C'
        complement['-']='-'
        complement['N']='N'
        revc=''
        for i in range(len(seq)):
            revc +=  complement[seq[i]]
        return(seq)
    
    def parseColorEncodingTable(self, color_file_name):
        with open(color_file_name, 'r') as color_file:
            for line in color_file:
                data = line.strip().split('\t')
                bp_2_seq = data[0].strip()
                pool = data[1].strip()
                spot = data[2].strip()
                #color = data[3].strip()
                code = data[4].strip()
                #color_to_num = self.colorToNumber(color)
                #tmp_key = pool + ';' + spot + ';' + color
                #tmp_key = pool + ';' + spot + ';' + color_to_num
                tmp_key = pool + ';' + spot + ';' + code
                if tmp_key in self.encoding_dic:
                    print "Duplicated color key", tmp_key
                self.encoding_dic[tmp_key] = bp_2_seq
        #print self.encoding_dic
    
    def parseCSVHeaders(self, fq_file, encoding_file_name):
        self.parseColorEncodingTable(encoding_file_name)
        for line in fq_file:
            if line.startswith('Feature'):
                data = line.strip().strip(',')
                data = data.split(',')
                self.idxToPoolList = [-1] * len(data)
                for i in xrange(len(data)):
                    if 'P' in data[i]:
                        pool = data[i].split('P')[-1].strip()
                        self.idxToPoolList[i] = pool
                return True
            else:
                print "Wrong CSV format"
                return False
        else: # No more lines to be read from file
            print "Empty CSV file"
            return False
            
    def decode_2mer(self, pool, spot, two_code):
        tmp_key = pool + ';' + str(spot) + ';' + two_code
        #print self.encoding_dic
        #print tmp_key
        if tmp_key in self.encoding_dic:
            #print '3333333', self.encoding_dic[tmp_key], tmp_key
            return self.encoding_dic[tmp_key]
        else:
            #print "No encoding information", tmp_key
            return 'NN'
        
    def decode(self, pool, code_str, molecule):
        if len(code_str) != 6:
            print "breakCode(): Wrong Code string", code_str
            return False
            
#        pool = '1'
#        code_str = '142413'
        
        molecule.cnt_reads += 1
        spot1 = code_str[0:2] # Spot1
        spot2 = code_str[2:4] # Spot2
        spot3 = code_str[4:6] # Spot3
        two_mer_1 = self.decode_2mer(pool, 1, spot1)
        two_mer_2 = self.decode_2mer(pool, 2, spot2)
        two_mer_3 = self.decode_2mer(pool, 3, spot3)
        cnt_x = 0
        if two_mer_1 == 'NN':
            cnt_x += 1
            two_mer_1 = 'NN'
        if two_mer_2 == 'NN':
            cnt_x += 1
            two_mer_2 = 'NN'
        if two_mer_3 == 'NN':
            cnt_x += 1
            two_mer_3 = 'NN'
        six_mer = two_mer_1 + two_mer_2 + two_mer_3
        six_mer_rev_comp = self.revcomp(six_mer)
#        print pool, code_str, spot1, spot2, spot3, six_mer, six_mer_rev_comp

                
        if cnt_x == 0:
            if param.pool_spec_orphans():
                if param.filterout_orphans() and six_mer_rev_comp not in self.all_sixmers[pool]:
                    if six_mer_rev_comp in self.nonsense_sixmers:
                        self.nonsense_sixmers[six_mer_rev_comp] += 1
                        self.cnt_filtered_out_sixmers += 1
                    else:
                        self.nonsense_sixmers[six_mer_rev_comp] = 1
                        self.cnt_filtered_out_sixmers += 1
                    return
                else:
                    #if param.filterout_orphans():
                    self.all_sixmers[pool][six_mer_rev_comp] += 1
                    self.cnt_used_sixmers += 1
                molecule.sixmer_list.append(six_mer_rev_comp)
            else:
                if param.filterout_orphans() and six_mer_rev_comp not in self.all_sixmers:
                    if six_mer_rev_comp in self.nonsense_sixmers:
                        self.nonsense_sixmers[six_mer_rev_comp] += 1
                        self.cnt_filtered_out_sixmers += 1
                    else:
                        self.nonsense_sixmers[six_mer_rev_comp] = 1
                        self.cnt_filtered_out_sixmers += 1
                    return
                else:
                    #if param.filterout_orphans():
                    self.all_sixmers[six_mer_rev_comp] += 1
                    self.cnt_used_sixmers += 1
                molecule.sixmer_list.append(six_mer_rev_comp)
        elif cnt_x == 1:
            molecule.sixmer_NN_list.append(six_mer_rev_comp)
        
    def parseCSVFile(self, fq_file, molecule_model_inst):
        for line in fq_file:
            if line.startswith('Features'):
                continue
            data = line.strip().split(',')
            feat = data[0].strip()
            x = data[1].strip()
            y = data[2].strip()
            pos_str = x + '.' + y + '_' + feat
            for i in xrange(3, len(data)):
                pool = self.idxToPoolList[i]
                if data[i] == '0': # No information
                    continue
                elif len(data[i]) == 4:
                    str_to_decode = '00' + data[i]
                elif len(data[i]) == 6:
                    str_to_decode = data[i]
                else: # ignore half colors and less than 2
                    continue
                self.decode(pool, str_to_decode, molecule_model_inst)
                #print len(molecule_model_inst.sixmer_list), len(molecule_model_inst.sixmer_NN_list)
            return True, pos_str
        else: # No more lines to be read from file
            return False, ''
        
    def parseBarcodeFile(self, fq_file, molecule_model_inst):
        for line in fq_file:
            data = line.strip().split('\t')
            gene_str = data[0].strip()
            sixmers_str = data[1].strip()
            sixmers_list = sixmers_str.split(',')
            pos_str = gene_str
            molecule_model_inst.sixmer_list = sixmers_list
            return True, pos_str
        else: # No more lines to be read from file
            return False, ''

    def checkInsufficientSixmers(self, molecule_model_inst):
        if len(molecule_model_inst.sixmer_list) < 1:
            del molecule_model_inst
            return True
        return False

    def checkInsufficientTargets(self, Assm_model_inst, molecule_model_inst, origin_gene):
        if origin_gene == '':
            #Assm_model_inst.num_molecules -= 1
            del molecule_model_inst
            return True
        Assm_model_inst.found_target_genes[origin_gene] += 1
        return False
        