# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 12:50:40 2017

@author: swoo
"""
import os
import Path_modules as Path
path = Path.Path()


# Parameters
ideal_cov_ = True # False # True
filterout_orphans_ = True # False # True
enable_known_mu_ = True # False # True
enable_blind_mu_ = True # False # True
show_XXX_targets_ = False # False # True

# Read masking
mask_ref_reads_ = True # False # True
mask_snp_reads_ = True # False # True
mask_mu_reads_ = True # False # True
mask_ref_N_thres_ = 1.0
mask_snp_N_thres_ = 1.0
mask_mu_N_thres_ = 1.0

mask_global_reads_ = True # False # True
mask_global_N_thres_ = 100.0

# 2-Spotters normalization
two_spotter_use_ = 1 # [0: dont use], [1 : divide 16], [2 : div locations]
dont_use_2_spot_in_global_ = False # False # True


# Using Genomic MTM score and garbage collectors for genomic dna data
is_genomic_dna_ = False # False # True


# MTM zscore
zscore_threshold_ = 2.5
zscore_delta_threshold_ = 0.5
print_raw_z_score_ = False # False # True


# Seperate ideal coverage plot
print_ideal_cov_ = True # False # True
ideal_file_name_ = os.path.join(path.encoding_dir, 'ideal_sixmers_1_fold.txt')


# Graph default coverage settings
default_ref_cov_ = 0.1
node_trimming_threshold_ = 0.5
default_mutation_cov_ = 0.0
fast_path_ = True # False # True


# USR Index File names
#mutation_file_name_ = 'predefined_mutations.txt'
#mutation_file_name_ = 'predefined_mutations_FFPE_run081.txt'
mutation_file_name_ = 'predefined_mutations_run083_artificial_100_test.txt'
USR_index_file_ = 'USR_index'
USR_mutation_index_file_ = 'USR_mutation_index'
USR_long_mutation_index_file_ = 'USR_long_mutation_index'
color_encoding_file_name_ = 'USR_color_table.txt'
mutation_file_name_ = os.path.join(path.input_dir, 'mutation', mutation_file_name_)
color_encoding_file_name_ = os.path.join(path.encoding_dir, color_encoding_file_name_)


# Output File names
assembled_read_out_ = '1_assembled_reads.fa'
consensus_graph_ = '2_consensus_graph_'
avg_consensus_graph_ = '3_avg_consensus_graph_'
max_consensus_graph_ = '5_max_consensus_graph_'
zscore_out_ = 'zscore.txt'
log_message_ = 'log.txt'
output_message_ = os.path.join(path.out_dir, 'output_log.txt')
debug_output_on_ = False # False # True
debug_output_file_name_ = os.path.join(path.out_dir, 'debug_output.txt')




# Parameter return functions
def output_message():
    return output_message_
def max_consensus_graph_out():
    return max_consensus_graph_
def dont_use_2_spot_in_global():
    return dont_use_2_spot_in_global_
def filterout_orphans():
    return filterout_orphans_    
def ideal_cov():
    return ideal_cov_    
def enable_known_mu():
    return enable_known_mu_    
def enable_blind_mu():
    return enable_blind_mu_    
def mask_ref_reads():
    return mask_ref_reads_    
def mask_snp_reads():
    return mask_snp_reads_
def mask_mu_reads():
    return mask_mu_reads_
def mask_ref_N_thres():
    return mask_ref_N_thres_ 
def mask_snp_N_thres():
    return mask_snp_N_thres_
def mask_mu_N_thres():
    return mask_mu_N_thres_
def two_spotter_use():
    return two_spotter_use_
def zscore_threshold():
    return zscore_threshold_
def zscore_delta_threshold():
    return zscore_delta_threshold_
def debug_output_on():
    return debug_output_on_
def get_debug_out_file_name():
    return debug_output_file_name_
def is_genomic_dna():
    return is_genomic_dna_    
def print_ideal_cov():
    return print_ideal_cov_
def print_raw_z_score():
    return print_raw_z_score_
def mutation_file_name():
    return mutation_file_name_
def color_encoding_file_name():
    return color_encoding_file_name_
def default_ref_cov():
    return default_ref_cov_
def trimming_threshold():
    return node_trimming_threshold_
def default_mutation_cov():
    return default_mutation_cov_
def ideal_file_name():
    return ideal_file_name_
def USR_index_file():
    return USR_index_file_
def USR_mutation_index_file():
    return USR_mutation_index_file_
def USR_long_mutation_index_file():
    return USR_long_mutation_index_file_
def consensus_graph_out():
    return consensus_graph_
def avg_consensus_graph_out():
    return avg_consensus_graph_
def zscore_out():
    return zscore_out_
def log_message():
    return log_message_
def assembled_read_out():
    return assembled_read_out_
def mask_global_reads():
    return mask_global_reads_
def mask_global_N_thres():
    return mask_global_N_thres_
def show_XXX_targets():
    return show_XXX_targets_
def fast_path():
    return fast_path_
    