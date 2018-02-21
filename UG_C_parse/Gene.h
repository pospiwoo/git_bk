//
// Created by swoo on 9/14/2017.
//
//#ifndef UGRAPHS_GENE_H
//#define UGRAPHS_GENE_H
//
//#endif //UGRAPHS_GENE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include "Graph.h"
#include "Tree.h"
#include "Parameters.h"
using namespace std;

struct ConstructGraph {
    Param* param;
    GGraph *GGraph_inst = new GGraph;
    string chr = "";
    string seq = "";
    string fasta_file;
    string test_file;
    string gff_file;
    string chr_in_process;

    ConstructGraph(){
        param = Param::getInstance();
        fasta_file = param->get_fasta_file();
        test_file = param->get_test_file();
        gff_file = param->get_gff_file();
        chr_in_process = param->get_chr_in_process();
    }

    void parseDNA(){
        string line;
        string header_str;
//        int cnt = 0;
        ifstream infile(fasta_file);
        while (getline(infile, line))
        {
            if (line[0] == '>')
                header_str = line;
            else
                seq = seq + line;
//            cnt++;
//            if(cnt%10000 == 1)
//                cout << cnt << ":" << seq.size() << endl;
        }
        infile.close();
    }

    void referenceG() {
        int from_idx = -1;
        for (int i = 0; i < seq.size(); i++) {
            GGraph_inst->addNode(i, seq[i], 0, 'R');
            GGraph_inst->addEdge(from_idx, i, 0);
            from_idx = i;
        }
    }

    void testReferenceG() {
        ofstream ofile(test_file);
        GraphNode *curr_node = GGraph_inst->NodeHash[-1];
        GraphNode *next_node = GGraph_inst->NodeHash[0];
        while(1){
            ofile << curr_node->base;
            if(curr_node->to_edges.size() == 0)
                break;
            else{
                next_node = curr_node->to_edges[0];
                curr_node = next_node;
            }
        }
        ofile.close();
    }

    void parseGFF(){
        ofstream ofile(test_file);
        string line;
        bool in_transcript = false;
        int node_from = -1;
        int tmp_cnt = 0;
        int tmp_line_cnt = 0;
        SuffixTree *tree = new SuffixTree();
        string trans_seq = "";
        string trans_id = "";
        string trans_id_cds = "";
        vector<int> list_coor_tree;
        vector<string> data;
        ifstream infile(gff_file);
        while (!infile.eof()){
            getline(infile, line);
            if (line[0] == '#' or line == "")
                continue;
            tokenizeLine(line, data);
            tmp_line_cnt++;

            string chr_str = data[0];
            if(chr_str.compare(chr_in_process) != 0)
                continue; // break
            string region_type = data[2];
            int start = stoi(data[3]) - 1;
            int end = stoi(data[4]);
            string tags = data[8];
//            cout << chr_in_process << ":" << start << "-" << end << " " << region_type << " ****** " << trans_id << endl;
            data.clear();

//            cout << tmp_cnt << endl;
//            cout << tmp_line_cnt << " : " << region_type << " " << start << " " << end << endl;

            if(region_type.compare("mRNA") == 0){
//                cout << list_coor_tree.size() << " : " << trans_seq << endl; /////////////////////////
                ofile  << trans_id << endl << trans_seq << endl;
                tokenizeTagInfo(tags, trans_id);
                in_transcript = true;
                node_from = -1;
                tmp_cnt += 1;
//                if(tmp_cnt > 1000)
//                    return;
                if(list_coor_tree.size() > 0) {
//                    cout << trans_seq << endl; /////////////////////////
                    tree->makeST(*tree, trans_seq, list_coor_tree);
                    list_coor_tree.clear();
                }
                trans_seq = "";
                continue;
            }
            else if(region_type.compare("exon") == 0){
                continue;
            }
            else if(region_type.find("UTR") != string::npos){
                continue;
            }
            else if(region_type.compare("CDS") == 0 && in_transcript == true){
                tokenizeTagInfoCDS(tags, trans_id_cds);
                if(trans_id.compare(trans_id_cds) != 0)
                    cout << "gff3 parsing is going wrong" << " " << trans_id_cds << " "<< trans_id << endl;
                string exon_seq = seq.substr(start,end-start); ////////////
                for(int j=0;j<exon_seq.size();j++){
                    int node_idx = start + j;
                    GGraph_inst->addNode(node_idx, exon_seq[j], 0, 'R');
                    GGraph_inst->addEdge(node_from, node_idx, 0);
                    node_from = node_idx;
                    list_coor_tree.push_back(node_idx);
                }
//                ofile << exon_seq << endl;
                trans_seq += exon_seq;
                continue;
            }
            else
                in_transcript = false;
        }
        infile.close();
        ofile.close();
    }

    void parseGFFtest(){
        cout << gff_file << endl<< endl;
        ofstream ofile(test_file);
        string line;
        bool in_transcript = false;
        int node_from = -1;
        int tmp_cnt = 0;
        string trans_seq = "";
        string trans_id = "";
        string trans_id_cds = "";
        vector<int> list_coor_tree;
        ifstream infile(gff_file);
        SuffixTree *tree = new SuffixTree();
        while (!infile.eof()){
            getline(infile, line);
            if (line[0] == '#' or line == "")
                continue;

            cout << "line:    " << line << endl;

            for(int i=0 ; i < line.size() ; i++ ){
                list_coor_tree.push_back(i);
            }
            tree->makeSTtest(*tree, line, list_coor_tree);
            list_coor_tree.clear();
            cout << endl << endl;
        }
        infile.close();
        ofile.close();
    }

    void tokenizeLine(string line, vector<string> &data) {
        istringstream  ss(line);
        std::string token;
        while(std::getline(ss, token, '\t'))
            data.push_back(token);
    }

    void tokenizeTagInfo(string tags, string &trans_id) {
        istringstream tagss(tags);
        string token1;
        vector<string> tags_vec;
        while (getline(tagss, token1, ';'))
            tags_vec.push_back(token1);
        for (int i = 0; i < tags_vec.size(); i++) {
            if (tags_vec[i].find("ID=transcript:") == 0) {
//                replace(tags_vec[i].begin(), tags_vec[i].end(), "ID=transcript:", "");
                trans_id = tags_vec[i].substr(14);
                return;
            }
//            else if (tags_vec[i].find("Parent=gene:") == 0) {
//                replace(tags_vec[i].begin(), tags_vec[i].end(), 'Parent=gene:', '');
//                string gene_id = tags_vec[i];
//            }
        }
        tags_vec.clear();
    }

    void tokenizeTagInfoCDS(string tags, string &trans_id_cds) {
        istringstream tagss(tags);
        string token1;
        vector<string> tags_vec;
        while (getline(tagss, token1, ';'))
            tags_vec.push_back(token1);
        for (int i = 0; i < tags_vec.size(); i++) {
            if (tags_vec[i].find("Parent=transcript:") == 0) {
//                replace(tags_vec[i].begin(), tags_vec[i].end(), 'Parent=gene:', '');
                trans_id_cds = tags_vec[i].substr(18);
//                cout << tags_vec[i] << " : " << trans_id_cds << endl;
                return;
            }
//            else if (tags_vec[i].find("Parent=gene:") == 0) {
//                replace(tags_vec[i].begin(), tags_vec[i].end(), 'Parent=gene:', '');
//                string gene_id = tags_vec[i];
//            }
        }
        tags_vec.clear();
    }

//        ofstream ofile(test_file);
//            for (int i=0; i< data.size(); i++)
//                ofile << "\t" << data[i];
//            ofile << endl;

//    ofile << ">" << header_str << endl;
//    ofile << all_str << endl;
//    for (std::map<int, GraphNode*>::iterator it=GGraph_inst->NodeHash.begin();
//         it!=GGraph_inst->NodeHash.end(); ++it)
//        delete it;
//    for(int i = 0; i < SNP.size(); i++) {
//        delete SNP[i];
//    }
//    for(int i = 0; i < SUB.size(); i++) {
//        delete SUB[i];
//    }
//    for(int i = 0; i < INS.size(); i++) {
//        delete INS[i];
//    }
//    for(int i = 0; i < DEL.size(); i++) {
//        delete DEL[i];
//    }
//    delete GGraph_inst;

};



