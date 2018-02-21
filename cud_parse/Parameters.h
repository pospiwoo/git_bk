//
// Created by swoo on 9/25/2017.
//

#ifndef DEV_C_0_01_PARAMETERS_H
#define DEV_C_0_01_PARAMETERS_H
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class Param {
    private:
        static Param* instance;
        bool enable_known_mu_ = false;
        bool enable_blind_mu_ = true;
        bool ideal_cov_ = false;
        bool filterout_orphans_ = false;
        bool pool_spec_orphans_ = false;
        bool show_XXX_targets_ = false;

        // Read masking
        bool mask_ref_reads_ = true;
        bool mask_snp_reads_ = true;
        bool mask_mu_reads_ = true;
        float mask_ref_N_thres_ = 1.0;
        float mask_snp_N_thres_ = 1.0;
        float mask_mu_N_thres_ = 1.0;

        bool mask_global_reads_ = true;
        float mask_global_N_thres_ = 100.0;

        // 2-Spotters normalization
        int two_spotter_use_ = 0; //# [0: dont use], [1 : divide 16], [2 : div locations]
        bool dont_use_2_spot_in_global_ = false;

        // Using Genomic MTM score and garbage collectors for genomic dna data
        bool is_genomic_dna_ = false;

        // MTM zscore
        float zscore_threshold_ = 2.5;
        float zscore_delta_threshold_ = 0.5;
        bool print_raw_z_score_ = false;

        // Seperate ideal coverage plot
        bool print_ideal_cov_ = false;
        //string ideal_file_name_ = os.path.join(path.encoding_dir, "ideal_sixmers_1_fold.txt")

        // Graph coverage thres settings
        float default_ref_cov_ = 0.1;
        float node_trimming_threshold_ = 0.5;
        float default_mutation_cov_ = 0.0;
        bool fast_path_ = true;

        // Input file names
//        string sequence_file_prefix = "C:\\Users\\swoo\\Box Sync\\dev\\precision_AGBT\\Run211";
//        string sequence_file_ = "Bioinformatics_S6_3Spot_C75_C174_Run211_cfDNA_HybSeq_1-100_SM75_174_filtered.csv";
//        string input_file_prefix = "C:\\Users\\swoo\\Box Sync\\dev\\Input\\";
//    string sequence_file_prefix = "/media/swoo/My Passport/Box Sync/dev/precision_AGBT/Run212/";
//    string sequence_file_ = "Bioinformatics_S6_3Spot_C61_C140_Run212_FFPE-4plex_HybSeq_61-140_SM61_140_filtered_combined.csv";
//    string input_file_prefix = "/media/swoo/My Passport/Box Sync/dev/precision_AGBT/encoding/";
//    string barcode_file_ = "barcode_list.txt"; // 'barcode_list_run190.txt' 'barcode_list.txt'
//    string target_file_ = "targets.txt"; // 'targets_KRAS.txt' 'targets_5'
//    string debug_file_ = "debug.txt";
        string sequence_file_prefix =
            "/home/swoo/Downloads/dev_cuda_0.01/Input/";
            //"/media/swoo/My Passport/dev_pyspark/Input/";
        string sequence_file_ =
            "med_runtime_100_freq_0_3.tsv";
            //"med_runtime_100_freq_0_3.tsv";
            //"med_runtime_10_freq_0_3_head_100.tsv";
        string input_file_prefix =
            "/home/swoo/Downloads/dev_cuda_0.01/Input/";
            //"/media/swoo/My Passport/dev_pyspark/Input/";
        string barcode_file_ = "barcode_list.txt"; // 'barcode_list_run190.txt' 'barcode_list.txt'
        string target_file_ = "target_sequences_tab_X.fa"; // 'targets_KRAS.txt' 'targets_5'
        string debug_file_ = "debug.txt";

        // USR Index File names
        string mutation_file_name_ = "predefined_mutations.txt";
        string USR_index_file_ = "USR_index";
        string USR_mutation_index_file_ = "USR_mutation_index";
        string USR_long_mutation_index_file_ = "USR_long_mutation_index";
        string color_encoding_file_name_ = "USR_color_table.txt";
        //string mutation_file_name_ = os.path.join(path.input_dir, "mutation", mutation_file_name_)
        //string color_encoding_file_name_ = os.path.join(path.encoding_dir, color_encoding_file_name_)

        // Output File names
        string assembled_read_out_ = "1_assembled_reads.fa";
        string consensus_graph_ = "2_consensus_graph_";
        string avg_consensus_graph_ = "3_avg_consensus_graph_";
        string max_consensus_graph_ = "5_max_consensus_graph_";
        string zscore_out_ = "zscore.txt";
        string log_message_ = "log.txt";
        //string output_message_ = os.path.join(path.out_dir, "output_log.txt")
        bool debug_output_on_ = false;
        //string debug_output_file_name_ = os.path.join(path.out_dir, "debug_output.txt")

        Param(); /* Private constructor to prevent instancing. */

    public:
        ofstream debugOutFile;
        static Param* getInstance(); /* Static access method. */
        bool is_enable_known_mu() { return enable_known_mu_; }
        bool is_enable_blind_mu() { return enable_blind_mu_; }
        bool is_ideal_cov() { return ideal_cov_; }
        bool is_filterout_orphans() { return filterout_orphans_; }
        bool is_pool_spec_orphans() { return pool_spec_orphans_; }
        bool is_show_XXX_targets() { return show_XXX_targets_; }
        bool is_mask_ref_reads() { return mask_ref_reads_; }
        bool is_mask_snp_reads() { return mask_snp_reads_; }
        bool is_mask_mu_reads() { return mask_mu_reads_; }
        float get_mask_ref_N_thres() { return mask_ref_N_thres_; }
        float get_mask_snp_N_thres() { return mask_snp_N_thres_; }
        float get_mask_mu_N_thres() { return mask_mu_N_thres_; }
        bool is_mask_global_reads() { return mask_global_reads_; }
        float get_mask_global_N_thres() { return mask_global_N_thres_; }
        int get_two_spotter_use() { return two_spotter_use_; }
        bool is_dont_use_2_spot_in_global() { return dont_use_2_spot_in_global_; }
        bool is_is_genomic_dna() { return is_genomic_dna_; }
        float get_zscore_threshold() { return zscore_threshold_; }
        float get_zscore_delta_threshold() { return zscore_delta_threshold_; }
        bool is_print_raw_z_score() { return print_raw_z_score_; }
        bool is_print_ideal_cov() { return print_ideal_cov_; }
        float get_default_ref_cov() { return default_ref_cov_; }
        float get_node_trimming_threshold() { return node_trimming_threshold_; }
        float get_default_mutation_cov() { return default_mutation_cov_; }
        bool is_fast_path() { return fast_path_; }
        string get_barcode_file() { return input_file_prefix + barcode_file_; }
        string get_target_file() { return input_file_prefix + target_file_; }
        string get_debug_file() { return input_file_prefix + debug_file_; }
        string get_sequence_file() { return sequence_file_prefix + sequence_file_; }
        string assembled_read_out() { return input_file_prefix + assembled_read_out_; }
    //    string get_() { return ; }
    //    string get_() { return ; }
    //    string get_() { return ; }
    //    string get_() { return ; }
    //
    //    // USR Index File names
    //    string mutation_file_name_ = "predefined_mutations.txt";
    //    string USR_index_file_ = "USR_index";
    //    string USR_mutation_index_file_ = "USR_mutation_index";
    //    string USR_long_mutation_index_file_ = "USR_long_mutation_index";
    //    string color_encoding_file_name_ = "USR_color_table.txt";
    //    string mutation_file_name_ = os.path.join(path.input_dir, "mutation", mutation_file_name_)
    //    string color_encoding_file_name_ = os.path.join(path.encoding_dir, color_encoding_file_name_)
    //
    //    // Output File names
    //    string assembled_read_out_ = "1_assembled_reads.fa";
    //    string consensus_graph_ = "2_consensus_graph_";
    //    string avg_consensus_graph_ = "3_avg_consensus_graph_";
    //    string max_consensus_graph_ = "5_max_consensus_graph_";
    //    string zscore_out_ = "zscore.txt";
    //    string log_message_ = "log.txt";
    //    string output_message_ = os.path.join(path.out_dir, "output_log.txt")
    //    bool debug_output_on_ = false;
    //    string debug_output_file_name_ = os.path.join(path.out_dir, "debug_output.txt")
    //
    //
    //    string get_chr_in_process() { return chr_in_process; }
    //    void put_chr_in_process(string chr_str) { chr_in_process = chr_str; }
    //    vector<string> get_all_chr() { return all_chr; }
};

Param* Param::instance = 0; /* Null, because instance will be initialized on demand. */
Param* Param::getInstance() {
    if (instance == 0) { instance = new Param(); }
    return instance;
}

Param::Param(){
    debugOutFile.open("/home/swoo/Downloads/debug.txt");
}

#endif //DEV_C_0_01_PARAMETERS_H
