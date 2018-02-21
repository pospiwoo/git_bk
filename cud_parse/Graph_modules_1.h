//
// Created by swoo on 9/14/2017.
//
#ifndef GRAPH_MODULES_H
#define GRAPH_MODULES_H
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <math.h>
//#include "Parameters.h"
using namespace std;

struct Node {
//    Param* param = Param::getInstance();
    int idx;
    char base;
    float cov;
    char type;
    bool visited = false;
    map<int, Node*> from_edges; // vector<GraphNode*>
    map<int, Node*> to_edges; // vector<GraphNode*>

    Node(int i, char b, int c, char t) {
        idx = i;
        base = b;
        cov = c;
        type = t;
    }
};



struct SpliceGraph {
    map<int, Node*> G; // map<int, GraphNode*> G;
    int graph_node_cnt;
    string max_seq = "";
    string max_seq_types = "";
    vector<float> max_seq_cov;
    vector<int> max_seq_score;
//    vector<int> DynamicP_track;
//    vector<float> DynamicP_score;

    SpliceGraph() {
        G[-1] = new Node(-1, 's', 0, 'S'); // G[0] = new GraphNode('s', 0, 'S');
        graph_node_cnt = 0;
    }
    ~SpliceGraph() {
        for (map<int, Node*>::iterator it = G.begin(); it != G.end(); ++it)
            delete it->second;
    }

    bool addNode(int idx, char base_str, float coverage, char type) {
        if (base_str == 'N') {
            cout << "this base is masked" << endl;
        }
        if (G.find(idx) == G.end()) {
            G[idx] = new Node(idx, base_str, coverage, type);
            graph_node_cnt += 1;
        }
        else {
            G[idx]->cov += coverage;
        }
        return true;
    }

    bool addEdge(int fro_idx, int to_idx, float coverage) {
        if (G.find(fro_idx) == G.end() ||
            G.find(to_idx) == G.end()) {
            cout << "Nodes not exist" << endl;
            return false;
        }
        if (G[fro_idx]->to_edges.find(to_idx) == G[fro_idx]->to_edges.end())
            G[fro_idx]->to_edges[to_idx] = G[to_idx]; // Assign edge from previous node
        if (G[to_idx]->from_edges.find(fro_idx) == G[to_idx]->from_edges.end())
            G[to_idx]->from_edges[fro_idx] = G[fro_idx]; // Assign edge to next node
        return true;
    }

    bool lookupMutationNode(int original_seq_len, int node_from, int node_to,
                char snp_char, int cov, char node_type) {
        if(graph_node_cnt == original_seq_len){ //Graph does not have any mutation node so we add
            int new_ind = graph_node_cnt;
            addNode(new_ind, snp_char, cov, node_type);
            addEdge(node_from, new_ind, cov);
            addEdge(new_ind, node_to, cov);
            return true;
        } else if(graph_node_cnt < original_seq_len){
            cout << "graph structure was not valid upon initiation" << endl;
            return false;
        }
        for (map<int, Node*>::iterator it = G.begin(); it != G.end(); ++it){
            int tmp_idx = it->first;
            Node *tmp_node = it->second;
            if(tmp_idx < original_seq_len)
                continue;
            if(tmp_node->from_edges.find(node_from) != tmp_node->from_edges.end()
               && tmp_node->to_edges.find(node_to) != tmp_node->to_edges.end()
                && tmp_node->base == snp_char)
                tmp_node->cov += cov;
                addEdge(node_from, tmp_idx, cov);
                addEdge(tmp_idx, node_to, cov);
                return true;
        }
        int new_ind = graph_node_cnt;
        addNode(new_ind, snp_char, cov, node_type);
        addEdge(node_from, new_ind, cov);
        addEdge(new_ind, node_to, cov);
        return true;
    }

    void grdQualityPath(int curr_idx) {
        Node *curr_node = G[curr_idx];
        curr_node->visited = 1;
        float max_cov = -1.f;
        int next_idx = -1;
        max_seq_types = max_seq_types + checkMutationString(curr_node);
        max_seq = max_seq + curr_node->base;
        max_seq_score.push_back(min(floor(curr_node->cov),9.0));
        max_seq_cov.push_back(curr_node->cov);
        if(curr_node->to_edges.size() > 0){
            for (map<int, Node*>::iterator it = curr_node->to_edges.begin();
                 it != curr_node->to_edges.end(); ++it){
                int to_idx = it->first;
                Node* next_node = it->second;
                if(next_node->visited == 0 && next_node->cov > max_cov){
                    max_cov = next_node->cov;
                    next_idx = to_idx;
                }
            }
            grdQualityPath(next_idx);
        }
        return;
    }

    char checkMutationString(Node *node){
        if(node->cov < 1.f){
            node->base = 'N';
            return 'u';
        }
        if(node->type == 'R')
            return '-';
        if(node->type == 'P')
            return '*';
        if(node->type == 'U')
            return '^';
        if(node->type == 'I')
            return '+';
        if(node->type == 'D')
            return '=';
        else
            return '?';
    }
    ///////////////////////////////////
//    bool addMutation(char base_str, int coverage, char type, int prev_ind,
//                     int next_ind, int ref_coor, int e_cov) {
//        for (int i = 0; i < NodeHash[prev_ind]->to_edges.size(); i++) {
//            if (NodeHash[prev_ind]->to_edges[i]->base == base_str &&
//                NodeHash[prev_ind]->to_edges[i]->type == type)
//                return false; // this mutation already exists
//        }
//        GraphNode *node_inst;
//        switch (type) {
//            case 'S': // SNP
//                node_inst = new GraphNode(SNP.size(), base_str, coverage, type);
//                SNP.push_back(node_inst);
//                break;
//            case 'U': // substitution
//                node_inst = new GraphNode(SUB.size(), base_str, coverage, type);
//                SUB.push_back(node_inst);
//                break;
//            case 'I': // insertion
//                node_inst = new GraphNode(INS.size(), base_str, coverage, type);
//                INS.push_back(node_inst);
//                break;
//            case 'D': // deletion
//                node_inst = new GraphNode(DEL.size(), base_str, coverage, type);
//                DEL.push_back(node_inst);
//                break;
//        }
//        addMutationEdge(node_inst, prev_ind, next_ind, e_cov);
//        graph_node_cnt += 1;
//        return true;
//    }
//
//    bool addMutationEdge(Node* node_inst, int fro_idx, int to_idx, int e_cov) {
//        NodeHash[fro_idx]->to_edges[to_idx] = NodeHash[to_idx]; // Assign edge from previous node
//        NodeHash[to_idx]->from_edges[fro_idx] = NodeHash[fro_idx]; // Assign edge to next node
//        return true;
//    }
//
//    pair< bool, vector<vector<int>> > iterateGenome(string search_seq, vector<int> coor) {
//        bool flag = false;
//        vector<int> single_list;
//        //        vector<vector<int>> return_list;
//        vector<vector<int> > return_list(10);
//        for (int i = 0; i<coor.size(); i++) {
//            //            cout << "qqweqweqwe111111221 : " << i << endl;
//            pair< bool, vector<int> > found_coor = retrieveGenome(search_seq, coor[i], 0);
//            //            cout << "11122222weqweqqweqweqwe111111221" << endl;
//            flag = found_coor.first;
//            single_list = found_coor.second;
//            return_list.push_back(single_list);
//        }
//        if (flag)
//            return make_pair(true, return_list);
//        else
//            return make_pair(false, return_list);
//    }
//
//    pair< bool, vector<int> > retrieveGenome(string search_seq, int coor_single, int error) {
//        vector<int> return_coor_list;
//        map<int, GraphNode*>::iterator it;
//
//        if (NodeHash.find(coor_single) == NodeHash.end()) {
//            //            cout << "11111111111" << endl;
//            return make_pair(false, return_coor_list);
//        }
//
//        GraphNode* cur_node = NodeHash[coor_single];
//
//        //        cout << "5151515 curr " << search_seq << " " << cur_node->base
//        //             << " " << cur_node->idx << endl;
//
//
//
//        if (search_seq.size() == 1 && cur_node->base == search_seq[0]) {
//            return_coor_list.push_back(cur_node->idx);
//            //            cout << "222222222 " << search_seq << " " << cur_node->base << " " << coor_single << endl;
//            return make_pair(true, return_coor_list);
//        }
//        else if (search_seq.size() == 1 && cur_node->base != search_seq[0]) {
//            //            cout << "21212122: " << cur_node->base << " " << search_seq << endl;
//            return make_pair(false, return_coor_list);
//        }
//            //        else if(cur_node->base == '')
//            //            return make_pair(false, return_coor_list);
//        else if (search_seq.size() > 1 && cur_node->to_edges.size() == 0) {
//            //            cout << "321131111231" << endl;
//            return make_pair(false, return_coor_list);
//        }
//
//
//        int i, j;
//        //        string tmp = "_";
//        //        for(i=0;i<cur_node->to_edges.size();i++){
//        //            tmp = tmp + cur_node->to_edges[i]->base + ":" + to_string(cur_node->to_edges[i]->idx) + "_";
//        //        }
//
//        for (it = cur_node->to_edges.begin(); it != cur_node->to_edges.end(); ++it) {
//            //            cout << "3333333 " << cur_node->to_edges.size()
//            //                    << " " << cur_node->idx
//            //                    << " " << it->second->idx << endl;
//
//            string next_str = search_seq.substr(1);
//            pair< bool, vector<int> > sub_result =
//                    retrieveGenome(next_str, it->second->idx, error);
//            bool flag = sub_result.first;
//            vector<int> sub_coor_list = sub_result.second;
//            //            if(flag == false)
//            //                continue;
//            if (flag == true) {
//                for (j = 0; j<sub_coor_list.size(); j++) {
//                    return_coor_list.push_back(sub_coor_list[j]);
//                }
//                return_coor_list.push_back(cur_node->idx);
//                return make_pair(true, return_coor_list);
//            }
//        }
//
//        if (return_coor_list.size() == 0)
//            return make_pair(false, return_coor_list);
//        else
//            return make_pair(true, return_coor_list);
//    }
};



#endif //GRAPH_MODULES_H