//
// Created by swoo on 9/14/2017.
//
//#ifndef UGRAPHS_GRAPH_H
//#define UGRAPHS_GRAPH_H
//
//#endif //UGRAPHS_GRAPH_H
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

struct GraphNode {
    int idx;
    char base;
    int cov;
    char type;
    bool visited = false;
    vector<GraphNode*> from_edges; // map<int, GraphNode*>
    vector<GraphNode*> to_edges; // map<int, GraphNode*>

    GraphNode(int i, char b, int c, char t){
        idx = i;
        base = b;
        cov = c;
        type = t;
    }
};


struct GGraph {
    map<int, GraphNode*> NodeHash; // map<int, GraphNode*> G;
    vector<GraphNode*> SNP;
    vector<GraphNode*> SUB;
    vector<GraphNode*> INS;
    vector<GraphNode*> DEL;
    int graph_node_cnt;

    GGraph(){
        NodeHash[-1] = new GraphNode(-1, 's', 0, 'S'); // G[0] = new GraphNode('s', 0, 'S');
        graph_node_cnt = 0;
    }

//    ~GGraph(){
//        for (std::map<int, GraphNode*>::iterator it=NodeHash.begin(); it!=NodeHash.end(); ++it)
//            delete it;
//        for(int i = 0; i < SNP.size(); i++) {
//            delete SNP[i];
//        }
//        for(int i = 0; i < SUB.size(); i++) {
//            delete SUB[i];
//        }
//        for(int i = 0; i < INS.size(); i++) {
//            delete INS[i];
//        }
//        for(int i = 0; i < DEL.size(); i++) {
//            delete DEL[i];
//        }
//    }

    bool addNode(int idx, char base_str, int coverage, char type){
        if(base_str == 'N'){
            cout << "this base is masked" << endl;
        }
        if(NodeHash.find(idx) == NodeHash.end()){
            NodeHash[idx] = new GraphNode(idx, base_str, coverage, type);
            graph_node_cnt += 1;
        }
        else{
            NodeHash[idx]->cov += coverage;
        }
        return true;
    }

    bool addEdge(int fro_idx, int to_idx, int coverage){
        if( NodeHash.find(fro_idx) != NodeHash.end() &&
            NodeHash.find(to_idx) != NodeHash.end() ){
            NodeHash[fro_idx]->to_edges.push_back(NodeHash[to_idx]); // Assign edge from previous node
            NodeHash[to_idx]->from_edges.push_back(NodeHash[fro_idx]); // Assign edge to next node
        }
        else{
            cout << "Nodes not exist" << endl;
            return false;
        }
        return true;
    }

    bool addMutation(char base_str, int coverage, char type, int prev_ind,
                     int next_ind, int ref_coor, int e_cov){
        for(int i = 0; i < NodeHash[prev_ind]->to_edges.size(); i++) {
            if(NodeHash[prev_ind]->to_edges[i]->base == base_str &&
               NodeHash[prev_ind]->to_edges[i]->type == type)
                return false; // this mutation already exists
        }
        GraphNode *node_inst;
        switch(type) {
            case 'S' : // SNP
                node_inst = new GraphNode(SNP.size(), base_str, coverage, type);
                SNP.push_back(node_inst);
                break;
            case 'U' : // substitution
                node_inst = new GraphNode(SUB.size(), base_str, coverage, type);
                SUB.push_back(node_inst);
                break;
            case 'I' : // insertion
                node_inst = new GraphNode(INS.size(), base_str, coverage, type);
                INS.push_back(node_inst);
                break;
            case 'D' : // deletion
                node_inst = new GraphNode(DEL.size(), base_str, coverage, type);
                DEL.push_back(node_inst);
                break;
        }
        addMutationEdge(node_inst, prev_ind, next_ind, e_cov);
        graph_node_cnt += 1;
        return true;
    }

    bool addMutationEdge(GraphNode* node_inst, int fro_idx, int to_idx, int e_cov){
        NodeHash[fro_idx]->to_edges.push_back(NodeHash[to_idx]); // Assign edge from previous node
        NodeHash[to_idx]->from_edges.push_back(NodeHash[fro_idx]); // Assign edge to next node
        return true;
    }
};



//#include <unordered_map>
//unordered_map<int, Node*> from_edges;
//unordered_map<int, Node*> to_edges;

////    struct Node node1;
//    Node* node1=new Node();
//    node1->idx = 100000;
//    node1->base = 'A';
//    node1->cov = 0;

////    struct Node node2;
//    Node* node2=new Node();
//    node2->idx = 200000;
//    node2->base = 'T';
//    node2->cov = 1;

//    struct Node node2;
//    strcpy( node1.subject, "C++ Programming");
//    Book1.book_id = 6495407;
//    cout << "IP Address: " << hashtable["www.element14.com"] << endl;
//    std::cout << "Hello, World!" << std::endl;
//    return 0;


//GraphNode *node1 = new GraphNode(100000, 'A', 0, 'R');
//GraphNode *node2 = new GraphNode(200000, 'T', 1, 'R');
//
////    node1->from_edges[0] = NULL;
//node1->to_edges[0] = node2;
//
//cout << "node1 to node2: " << node1->idx << " -> " << node1->base << endl;

//        vector<GraphNode*>::iterator v = NodeHash[prev_ind]->to_edges.begin();
//        while( v != NodeHash[prev_ind]->to_edges.end()) {}//
// Created by swoo on 9/14/2017.
//
