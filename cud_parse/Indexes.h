//
// Created by swoo on 9/26/2017.
//
#ifndef DEV_C_0_01_INDEXES_H
#define DEV_C_0_01_INDEXES_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <algorithm>
#include <cctype>
#include <locale>
#include "Parameters.h"
using namespace std;

class Indexes{
private:
    Param* param = Param::getInstance();
    static Indexes* instance;
    string target_file;

    Indexes(); /* Private constructor to prevent instancing. */

public:
    map<string, map<string, vector<int> > > sixmer_dic_target;
    map<string, string> target_sequences;
    map<string, float> MTM_score;
    map<char, char> comp_map;

    static Indexes* getInstance(); /* Static access method. */

    Param* get_param() {
        return param;
    }

    void buildIndexes(){
        target_file = param->get_target_file();
        string line;
        vector<string> data;
        // >BRAF_1
        // TTGTAGACTGTTCCAAATGATCCAGATCCAATTCTTTGTCCCACTGTAATCTGCCCATCAGGAATCT...
        ifstream inTargetFile(target_file.c_str());
        while (!inTargetFile.eof()){
            getline(inTargetFile, line);
            if (line[0] == '#' || line == "")
                continue;
            tokenizeLine(line, data);

            string gene_str = data[0];
            string seq = data[3];
            transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
            data.clear();
//            cout << line << endl << gene_str << endl << seq << endl << endl;
            target_sequences[gene_str] = seq;
            MTM_score[gene_str] = 0.0f;
            for(int i=0; i<seq.size()-5; i++){
                string key_seq = seq.substr(i,6);
                if (sixmer_dic_target.find(key_seq) == sixmer_dic_target.end()) {
                    vector<int> tmp_idx_list;
                    tmp_idx_list.push_back(i);
                    map<string, vector<int> > tmp_dic;
                    tmp_dic[gene_str] = tmp_idx_list;
                    sixmer_dic_target[key_seq] = tmp_dic;
                } else{
                    if (sixmer_dic_target[key_seq].find(gene_str) == sixmer_dic_target[key_seq].end()){
                        sixmer_dic_target[key_seq][gene_str].push_back(i);
                    } else{
                        vector<int> tmp_idx_list;
                        tmp_idx_list.push_back(i);
                        sixmer_dic_target[key_seq][gene_str] = tmp_idx_list;
                    }
                }
            }
        }
//        print_target_index();
    }

    void tokenizeLine(string line, vector<string> &data) {
        istringstream  ss(line);
        string token;
        while(getline(ss, token, '\t'))
            data.push_back(token);
    }

    void tokenizeLine(string line, vector<string> &data, char delimiter) {
        istringstream  ss(line);
        string token;
        while(getline(ss, token, delimiter))
            data.push_back(token);
    }

    string revComp(string seq) {
        string rev_comp = "";
        if(seq.size()==1)
            rev_comp = comp_map[seq[0]];
        else
            for(int i=0; i<seq.size(); i++){
                rev_comp =  comp_map[seq[i]] + rev_comp;
            }
        return rev_comp;
    }

    string convert_seq(string seq) {
        string new_seq = "";
        for(int i=0; i<seq.size(); i++){
            if(seq[i] == '2')
                new_seq += "G";
            else if(seq[i] == '4')
                new_seq += "R";
            else if(seq[i] == '3')
                new_seq += "Y";
            else if(seq[i] == '1')
                new_seq += "B";
        }
        return new_seq;
    }
    
    /*
    static inline void ltrim(string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) {
            return !isspace(ch);
        }));
    }

    static inline void rtrim(string &s) { // trim from end (in place)
        s.erase(find_if(s.rbegin(), s.rend(), [](int ch) {
            return !isspace(ch);
        }).base(), s.end());
    }

    static inline void trim(string &s) { // trim from both ends (in place)
        ltrim(s);
        rtrim(s);
    }

    static inline string ltrim_copy(string s) { // trim from start (copying)
        ltrim(s);
        return s;
    }

    static inline string rtrim_copy(string s) { // trim from end (copying)
        rtrim(s);
        return s;
    }

    static inline string trim_copy(string s) { // trim from both ends (copying)
        trim(s);
        return s;
    }
*/

    void print_target_index(){
        map<string, map<string, vector<int> > >::iterator it;
        map<string, vector<int> >::iterator it_1;

        for (it = sixmer_dic_target.begin(); it != sixmer_dic_target.end(); ++it) {
            string key_seq = it->first;
            map<string, vector<int> > gene_data = it->second;
            cout << key_seq << " - " ;
            for (it_1 = gene_data.begin(); it_1 != gene_data.end(); ++it_1) {
                string gene_str = it_1->first;
                vector<int> idx_list = it_1->second;
                cout << gene_str << ":" ;
                if(idx_list.size() == 1)
                    cout << idx_list[0];
                else
                    for(int i=0; i<idx_list.size(); i++)
                        cout << idx_list[i] << ",";
                cout << " ";
            }
            cout << endl;
        }
    }

//    map<string, string> get_default_ref_cov() { return target_sequences; }
};

Indexes* Indexes::instance = 0; /* Null, because instance will be initialized on demand. */
Indexes* Indexes::getInstance() {
    if (instance == 0) { instance = new Indexes(); }
    return instance;
}

Indexes::Indexes(){
    comp_map['A'] = 'T';
    comp_map['T'] = 'A';
    comp_map['C'] = 'G';
    comp_map['G'] = 'C';
    comp_map['-'] = '-';
    comp_map['N'] = 'N';
}

#endif //DEV_C_0_01_INDEXES_H
