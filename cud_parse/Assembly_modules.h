//
// Created by swoo on 9/25/2017.
//

#ifndef DEV_C_0_01_ASSEMBLY_MODULES_H
#define DEV_C_0_01_ASSEMBLY_MODULES_H
#include "Graph_modules.h"
#include "Parameters.h"
#include "Indexes.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
using namespace std;

struct MoleculeModel{
    Param* param = Param::getInstance();
    vector<string> sixmer_list;
//    vector<string> sixmer_N_list;
//    vector<string> sixmer_NN_list;
    vector<string> mutation_search_list;
//    vector<string> mutation_barcode_list;
    SpliceGraph* Graph;
    Indexes *index;
    string gene_str = "";
    string pos_header = "";
    float default_ref_cov;
    int cnt_reads = 0;
    int gene_uniq_sixmer = 0;
    int mutation_barcode_cnt = 0;

    MoleculeModel(string gene_str_v){
        index = Indexes::getInstance();
        default_ref_cov = param->get_default_ref_cov();
        gene_str = gene_str_v;
        initGraph();
    }
    ~MoleculeModel(){
        delete Graph;
    }

    void initGraph(){
        Graph = new SpliceGraph();
        string seq_to_search = index->target_sequences[gene_str];
        for(int i=0; i<seq_to_search.size(); i++){
            Graph->addNode(i, seq_to_search[i], default_ref_cov, 'R');
        }
        for(int i=1; i<seq_to_search.size()-1; i++){ //for i in xrange(1, len(seq_to_search)-1){
            Graph->addEdge(i-1, i, default_ref_cov);
            Graph->addEdge(i, i+1, default_ref_cov);
        }
    }

//    void initGraphCopy(){
//        Graph = index->copyGraph(Graph, gene_str);
//    }

    void addSNP(int position, char snp_char, int cov){
        int from_node = position - 1;
        int to_node = position + 1;
        int original_seq_len = index->target_sequences[gene_str].size();
        Graph->lookupMutationNode(original_seq_len, from_node, to_node,
                                 snp_char, cov, 'P');
    }
};



struct Assembly{
    Param* param;
    Indexes *index;
//    Align_module = Align.Alignment()
    map<string, string> barcode_dic;
    map<string, bool> all_barcode_list;
    map<int, tuple<char, string>> snp_dic;
    vector<string> cur_sixmer_list;
    vector<string> cur_mut_sixmer_list;
    float default_ref_cov;
    float trimming_threshold;
    int num_molecules;
    int num_reads;
    string sep;
    string xy_str;

    Assembly(){
        param = Param::getInstance();
        index = Indexes::getInstance();
//    Align_module = Align.Alignment()
        default_ref_cov = param->get_default_ref_cov();
        trimming_threshold = param->get_node_trimming_threshold();
        num_molecules = 0;
        num_reads = 0;
        sep = "\t";
//        parseBarcode();
    }
    ~Assembly(){
        param->debugOutFile.close();
    }

    void parseBarcode() {
        //TGCGTG	RRGGBR
        //GCTTGA	BYBGYR
        string barcode_file = param->get_barcode_file();
        vector<string> data;
        string line;
        ifstream inBarcodeFile(barcode_file);
        while (!inBarcodeFile.eof()) {
            getline(inBarcodeFile, line);
            if (line[0] == '#' || line == "")
                continue;
            index->tokenizeLine(line, data);
            string raw_seq = index->rtrim_copy(data[0]);
            string seq = index->revComp(raw_seq);
            string color = index->rtrim_copy(data[1]);
            data.clear();
            barcode_dic[color] = seq;
            all_barcode_list[seq] = true;
        }
    }

    void process() {
        string sequence_file = param->get_sequence_file();
        ifstream inSequnceFile(sequence_file);
        bool flag_keep_reading = true;
        string origin_gene;
        while(flag_keep_reading == true){
//            flag_keep_reading = parseCSVFile(inSequnceFile);
            flag_keep_reading = parseTSVFile(inSequnceFile);
            origin_gene = mapToTargets();
            if(origin_gene == "")
                continue;
            MoleculeModel *molecule_model = new MoleculeModel(origin_gene);
//            molecule_model.initMutationGraph(origin_gene)
            estimateAllCov(molecule_model, origin_gene);
//            self.updateReadCounts(molecule_model, xy_str)
            FindMutations(molecule_model, origin_gene);
            molecule_model->Graph->grdQualityPath(0);
            PrintFullReads(molecule_model, "", origin_gene);

//            for (map<int, Node*>::iterator itg = molecule_model->Graph->G.begin();
//                 itg != molecule_model->Graph->G.end(); ++itg){
//                int tmp_int = itg->first;
//                Node *tmp_node = itg->second;
//                param->debugOutFile << floor(tmp_node->cov) << " ";
//            }
//            param->debugOutFile << endl;

            clearMEM();
            delete molecule_model;
        }
    }

    bool PrintFullReads(MoleculeModel *molecule_inst, string header_str, string gene_str) {
        param->debugOutFile << ">" << gene_str << ";" << header_str << endl;
        param->debugOutFile << molecule_inst->Graph->max_seq << endl;
        param->debugOutFile << molecule_inst->Graph->max_seq_types << endl;
        for(int i=0; i<molecule_inst->Graph->max_seq_score.size(); i++){
            param->debugOutFile << molecule_inst->Graph->max_seq_score[i];
        }
        param->debugOutFile << endl;
    }

    void FindMutations(MoleculeModel *molecule_inst, string gene_str){
        if(cur_mut_sixmer_list.size() == 0)
            return;
        for(int i=0; i<cur_mut_sixmer_list.size(); i++){
            ProcessHamming(molecule_inst,
                      index->target_sequences[gene_str],
                      cur_mut_sixmer_list[i]);
            if(snp_dic.size() > 0)
                iterSNP(molecule_inst, gene_str, index->target_sequences[gene_str].size());
        }
    }

    void ProcessHamming(MoleculeModel *molecule_inst, string target_seq, string sixmer){
        int i;
        for(i=0; i<target_seq.size()-6; i++){
            string target_chunk = target_seq.substr(i,6);
            hammingDist(target_chunk, sixmer);
        }
    }

    void hammingDist(string target_str, string six_str){
        int diffs = 0;
        int last_diff_idx;
        int i;
        char ref_char, diff_char;
        if(target_str.size() == six_str.size()){
            for(i=0; i<target_str.size(); i++){
                if(target_str[i] != six_str[i]){
                    diffs++;
                    last_diff_idx = i;
                    ref_char = target_str[i];
                    diff_char = six_str[i];
                }
            }
        }
        else
            cout << "hammingDist(): length different 1 " << target_str << " 2 " << six_str << endl;
        if(diffs == 1){
            tuple<char, string> tuple_tmp = make_tuple(diff_char, target_str);
            snp_dic[last_diff_idx] = tuple_tmp;
//            param->debugOutFile << diffs << " "
//                                << last_diff_idx << " "
//                                << ref_char << "->" <<
//                                diff_char << ":" << target_str << " "
//                                << six_str << endl;
        }
    }

    void iterSNP(MoleculeModel *molecule_inst, string gene_str, int target_seq_len){
        int i, tmp_pos, snp_idx;
        float normalized_vote = 1.f / (float)snp_dic.size();
        map<int, tuple<char, string>>::iterator it;
        for (it = snp_dic.begin(); it != snp_dic.end(); ++it){
            int rel_snp_idx = it->first;
            tuple<char, string> tuple_tmp = it->second;
            char snp_char = get<0>(tuple_tmp);
            string target_chnk = get<1>(tuple_tmp);
            vector<int> locations = index->sixmer_dic_target[target_chnk][gene_str];
            for(i=0; i<locations.size(); i++){
                tmp_pos = locations[i];
                snp_idx = tmp_pos + rel_snp_idx;
                if(snp_idx <= 2 || snp_idx >= (target_seq_len-3))
                    continue;
            }
            molecule_inst->addSNP(snp_idx, snp_char, normalized_vote);
            //self.Assm_model.addGlobalSNP(gene_str, snp_idx, snp_char, normalized_vote)
            increaseCov(molecule_inst, gene_str, tmp_pos, normalized_vote, rel_snp_idx);
        }
    }

    void estimateAllCov(MoleculeModel *molecule_model, string gene_str){
        Coverage(molecule_model, gene_str); // Perfect match
//        self.Assm_model.MutationModel.Coverage(molecule_model, origin_gene) // Predefined short mutations
//        self.Assm_model.MutationModel.CoverageLong(molecule_model, origin_gene) // Predefined Long mutations
    }

    void Coverage(MoleculeModel *molecule_model, string gene_str){
        for(int i=0; i<cur_sixmer_list.size(); i++){
            string tmp_sixmer = cur_sixmer_list[i];
            if(index->sixmer_dic_target.find(tmp_sixmer) !=
               index->sixmer_dic_target.end()
               && index->sixmer_dic_target[tmp_sixmer].find(gene_str) !=
               index->sixmer_dic_target[tmp_sixmer].end()){
                float normalized_vote = 1.f / (float)index->sixmer_dic_target[tmp_sixmer][gene_str].size();
                iterLocation(molecule_model, gene_str, tmp_sixmer, normalized_vote);
            } else
                cur_mut_sixmer_list.push_back(tmp_sixmer);
        }
    }

    void iterLocation(MoleculeModel *molecule_model, string gene_str, string sixmer, float normalize){
        int i = 0, tmp_pos = -1;
        vector<int> tmp_vec = index->sixmer_dic_target[sixmer][gene_str];

//        index->print_target_index(); ///////////////////////////////////////
//        cout << index->sixmer_dic_target["TTTTGG"]["APC_3"][0];

//        param->debugOutFile << "1111111" << endl;
//        param->debugOutFile << index->sixmer_dic_target[sixmer][gene_str][0] << endl;
//        param->debugOutFile << "1111111" << endl;
//
//        param->debugOutFile << gene_str << " " << sixmer << " "
//                << tmp_vec.size() << " "
//                << tmp_vec[0]
//                << ": ";
        for(i=0; i < tmp_vec.size(); i++){
//            param->debugOutFile << i << "-" << tmp_vec[i] << " ";
            increaseCov(molecule_model, gene_str, tmp_vec[i], normalize);
        }
//        param->debugOutFile << endl;
//        cout << index->sixmer_dic_target[sixmer][gene_str].size() << endl;
    }

    void increaseCov(MoleculeModel *molecule_model, string gene_str, int start_idx, float normalize){
        for(int i=0; i<6; i++){
            molecule_model->Graph->G[start_idx+i]->cov += normalize;
//            molecule_model.g_Graph[gene_str].G[start_idx+i].cov += normalized_vote // Update Global Graph
        }
    }

    void increaseCov(MoleculeModel *molecule_model, string gene_str, int start_idx, float normalize, int dont_increase){
        for(int i=0; i<6; i++){
            if(i != dont_increase)
                molecule_model->Graph->G[start_idx+i]->cov += normalize;
        }
    }

    string mapToTargets(){
        string match_gene_str = "";
        map<string, vector<int>>::iterator it_1;
        for(int i=0; i<cur_sixmer_list.size(); i++){
            string tmp_sixmer = cur_sixmer_list[i];
            if(index->sixmer_dic_target.find(tmp_sixmer)
               != index->sixmer_dic_target.end()){
                map<string, vector<int>> gene_data = index->sixmer_dic_target[tmp_sixmer];
                for (it_1 = gene_data.begin(); it_1 != gene_data.end(); ++it_1){
                    string gene_str = it_1->first;
                    index->MTM_score[gene_str] += 1.0;
                }
            }
        }

        map<string, float>::iterator it;
        float max_score = 0.0f;
        for (it = index->MTM_score.begin(); it != index->MTM_score.end(); ++it){
//            param->debugOutFile << it->first << ":" << it->second << " ";
            if(it->second > max_score){
                max_score = it->second;
                match_gene_str = it->first;
            }
        }
        if(max_score == 0.0f)
            return "";
        return match_gene_str;
    }

    void clearMEM(){
        cur_sixmer_list.clear();
        cur_mut_sixmer_list.clear();
        clearMTMScore();
        clear_snp_dic();
    }

    void clearMTMScore(){
        map<string, float>::iterator it;
        for (it = index->MTM_score.begin(); it != index->MTM_score.end(); ++it) {
            it->second = 0.0f;
        }
    }

    void clear_snp_dic(){
        map<int, tuple<char, string>>::iterator it;
        for (it = snp_dic.begin(); it != snp_dic.end(); ++it) {
            snp_dic.erase(it);
        }
    }

    bool parseCSVFile(ifstream& inSequnceFile) {
        vector<string> data;
        string line;
        int i, cnt_zero = 0;
        if(inSequnceFile.eof())
            return false;
        getline(inSequnceFile, line);
        if (line[0] == '#' || line == "")
            return false;
        if (line[0] == 'F') // line.startswith('Features'):
            return true;
        index->tokenizeLine(line, data, ',');
        if (data.size() < 5)
            return false;
        string features = data[0];
        string fov = data[1];
        string x = data[2];
        string y = data[3];
        xy_str = fov + "_" + x + "." + y;
//        param->debugOutFile << xy_str << "=" << data.size() << " : ";
        for (i=4; i<data.size(); i++) {
            if (data[i] == "0" || data[i].size() != 6){
                cnt_zero++;
                continue;
            }
            string conseq1 = index->convert_seq(data[i]);
            string seq = barcode_dic[conseq1];
            if(seq.size() == 6)
                cur_sixmer_list.push_back(seq);
//            param->debugOutFile << conseq1 << ";" << seq << " ";
        }
//        param->debugOutFile << data.size() << "-" << cnt_zero << " = ";
//            if param.remove_sticky_barcodes():
//                sixmer_list_tmp = self.remove_sticky(sixmer_list_tmp)
        data.clear();
        return true;
    }

    bool parseTSVFile(ifstream& inSequnceFile) {
        vector<string> data;
        string line;
        int i, cnt_zero = 0;
        if(inSequnceFile.eof())
            return false;
        getline(inSequnceFile, line);
        if (line[0] == '#' || line == "")
            return false;
        if (line[0] == 'F') // line.startswith('Features'):
            return true;
        index->tokenizeLine(line, data, '\t');
        if (data.size() < 5)
            return false;
        string features = data[0];
        string fov = data[1];
        string x = data[2];
        string y = data[3];
        xy_str = fov + "_" + x + "." + y;
//        param->debugOutFile << xy_str << "=" << data.size() << " : ";
        for (i=4; i<data.size(); i++) {
            if (data[i] == "0" || data[i].size() != 6){
                cnt_zero++;
                continue;
            }
            string seq = data[i];
            if(seq.find('N')!=string::npos)
                continue;
            if(seq.size() == 6)
                cur_sixmer_list.push_back(seq);
//            param->debugOutFile << conseq1 << ";" << seq << " ";
        }
//        param->debugOutFile << data.size() << "-" << cnt_zero << " = ";
//            if param.remove_sticky_barcodes():
//                sixmer_list_tmp = self.remove_sticky(sixmer_list_tmp)
        data.clear();
        return true;
    }
};


#endif //DEV_C_0_01_ASSEMBLY_MODULES_H


//    string mutation_file_name = param.mutation_file_name()
//    color_encoding_file_name = param.color_encoding_file_name()
//    fastq_count_name = fastq_count
//    fastq_name = fastq
//    target_fa_name = target_fa
//    target_sequences = {}
//    g_Graph = {}
//    found_target_genes = defaultdict(int)
//    buildUSRindex()
//    Graph_view = Graph.GraphView(Assm_view.path)
//    MutationModel = Mutation.MutationModel(self, Graph_view)
//    buildGlobalMutationGraph()  // Initiate global Mutation Graph
//    ideal_barcode_cov = {}

//    vector<string> NN_subst = ["AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG","CA", "CT", "CC", "CG","GA", "GT", "GC", "GG"];
//    initMiscVars();
//    int cnt_sixmer = 0;
//    int cnt_no_match_sixmer = 0;
//    int cnt_NN_mer = 0;
//    sum_zscores = []
//    cnt_zscores = 0
//    delta_zscore = 0.0
//    MTM_removed_cnt = 0
//    XXX_removed_cnt = 0
//    MTM_score = {}

//    for gene_str in Assm_model.target_sequences:
//            MTM_score[gene_str] = 0.0
//    sum_zscores.append(0.0)
//    if param.debug_output_on():
//            debugFile = open(param.get_debug_out_file_name(),'w')