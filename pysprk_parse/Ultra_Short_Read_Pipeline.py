# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 09:08:20 2016

@author: Sunghee Woo
"""
import sys, os, time
import Assembly_modules as Assm
import Path_modules as Path

class USRController:
    def __init__(self, fastq, target_fa, path_inst):
        self.Assm_view = Assm.AssemblyView(path_inst)
        self.Assm_model = Assm.AssemblyModel(fastq, target_fa, self.Assm_view)
        self.Assm_controller = Assm.AssemblyController(self.Assm_model, self.Assm_view)
        
    def Process(self):
        start = time.time()
        self.Assm_controller.Process()
        end = time.time()
        self.Assm_view.print_time(start, end)
        

if __name__=='__main__':
    path = Path.Path()
    if len(sys.argv) == 1:
        
        fastq = 'sample_2_low.fastq'
        fastq = 'sample_2_med.fastq'
        fastq = 'sample_2_high.fastq'
        
        fastq = 'sample_2_low.fastq'
        target_fa = 'target_sequeces.fa'
        
        target_fa = 'target_sequeces_mod.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_subs.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_subs_long.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_del_long.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_del.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_ins_long.fa'
        target_fa = 'target_sequeces_mod_BRAF_1_ins.fa'
        
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_F1_old.csv'
        target_fa = 'target_sequences_barcode_test_single.fa'
        
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_all.csv'
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_F1.csv'
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_F1_20161219_half_cycle.csv'
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_all_20161219.csv'
        
        target_fa = 'target_sequences_barcode_test_BRAF.fa'
        target_fa = 'target_sequences_barcode_test_arti_mod.fa'
        target_fa = 'target_sequences_barcode_test_arti_mod_2.fa'
        target_fa = 'target_sequences_barcode_test_arti_mod_3.fa'
        target_fa = 'target_sequences_barcode_test_with_SNPs.fa'
        target_fa = 'target_sequences_barcode_test_arti_mod.fa'
        
        fastq = 'Run_old\\Bioinformatics_S6_SEQ_HighConf_nH2_F1_0.csv'
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_all_20161228.csv'
        fastq = 'Bioinformatics_S6_SEQ_HighConf_nH2_Run073_20161230.csv'
        fastq = 'Run_old\\Bioinformatics_S6_SEQ_HighConf_nH2_D2_F3.csv'
        fastq = 'Run076\\Bioinformatics_S6_SEQ_HighConf_nH2_F4.csv'
        fastq = 'Run075\\Bioinformatics_S6_SEQ_HighConf_nH2_F3.csv'
        fastq = 'Run074\\Bioinformatics_S6_SEQ_HighConf_nH2_F3.csv'
        fastq = 'Run076_8_pools\\Bioinformatics_S6_SEQ_HighConf_nH5_F5.csv'
        fastq = 'Run076_8_pools\\Bioinformatics_S6_SEQ_HighConf_nH5_all_FOV.csv'
        fastq = 'Run_old\\Bioinformatics_S6_SEQ_HighConf_nH2_D2_all_20161228.csv'
        fastq = 'Run074_100_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_all.csv'
        fastq = 'Run079\\Bioinformatics_S6_SEQ_HighConf_nH5_all.csv'
        fastq = 'Run080\\SuperMask_Data_5min_Hyb\\Bioinformatics_S6_SEQ_HighConf_nH2_all.csv'
        fastq = 'Run080\\SuperMask_Data_1min_Hyb\\Bioinformatics_S6_SEQ_HighConf_nH2_all.csv'
        fastq = 'Run081\SuperMask_Data\\Bioinformatics_S6_SEQ_HighConf_nH5_all.csv'
        fastq = 'Run081\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_all.csv'
        fastq = 'Run075_100_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_F3.csv'
        fastq = 'Run082_50_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_all.csv'
        fastq = 'Run075_100_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_all.csv'
        fastq = 'Run079\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_100_cycle.csv'
        fastq = 'Run083_50_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_all.csv'
        fastq = 'Run082_100_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_all_D2.csv'
        fastq = 'Run083_75_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_all_D2.csv'
        fastq = 'Run085\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_all.csv'
        fastq = 'Run083_75_cycle\\Bioinformatics_S6_SEQ_HighConf_nH5_all_D2_p123.csv'
        fastq = 'Run085\\Bioinformatics_S6_SEQ_HighConf_nH5_D2_all_p123.csv'
        
        target_fa = 'target_sequences_genomic_dna.fa'
        target_fa = 'target_sequences_genomic_dna_XXX.fa'
        target_fa = 'target_sequences_barcode_test_short.fa'
        target_fa = 'target_sequences_genomic_dna_XXX_sep_wrong.fa'
        #target_fa = 'target_sequences_barcode_KRAS_COSM518_diff.fa'
        #target_fa = 'target_sequences_barcode_KRAS_COSM522_diff.fa'
        #target_fa = 'target_sequences_barcode_KRAS_COSM532_diff.fa'
        target_fa = 'target_sequences_genomic_dna_XXX_sep.fa'
        target_fa = 'target_sequences_barcode_test_artificial_100_mutation.fa'
        target_fa = 'target_sequences_barcode_test.fa'
        
    
    elif len(sys.argv) < 3:
        print('Insufficient arguments')
        print('Usage: python XXX.py <FASTQ> <Target_FASTA>')
        exit()
    elif len(sys.argv) >= 3:
        fastq = sys.argv[1]
        target_fa = sys.argv[2]
    
    
    fastq = os.path.join(path.input_dir, fastq)
    target_fa = os.path.join(path.input_dir, target_fa)
    
    Assm_view = Assm.AssemblyView(path)
    Assm_model = Assm.AssemblyModel(fastq, target_fa, Assm_view)
    Assm_controller = Assm.AssemblyController(Assm_model, Assm_view)
    
    start = time.time()
    Assm_controller.Process() # Core module
    end = time.time()
    Assm_view.print_time(start, end)
