//
// Created by swoo on 9/15/2017.
//
//#ifndef UGRAPHS_TREE_H
//#define UGRAPHS_TREE_H
//
//#endif //UGRAPHS_TREE_H
#include <iostream>
#include <string>
#include <vector>
#include <map>
//#include "Parameters.h"
using namespace std;

struct TreeNode {
    string lab;
    map<int, bool> coor;
    map<char, TreeNode*> out; // vector<TreeNode*> out;

    TreeNode(string seq_str, int i){
        lab = seq_str;
        coor[i] = true;
    }
};
struct SuffixTree {
//    Param* param;
    TreeNode *root;
    int tree_node_len;
    int MAX_IDX;
    int MIN_IDX;

    SuffixTree(){
        root = new TreeNode("s", -1);
        tree_node_len = 0;
//        param = Param::getInstance();
//        MAX_IDX = param->get_MAX_IDX();
        MAX_IDX = 110;
        MIN_IDX = 10;
    }

    bool makeST(SuffixTree &tree, string seq_str, vector<int> coor_list){
        if (seq_str.size() < MIN_IDX || seq_str == "")
            return false;
        if(seq_str.size() != coor_list.size()){
            cout << "seq and coor_list length is different" << endl;
            return false;
        }
        int i;
        for(i = 0; i < seq_str.size() - MAX_IDX ; i++){
            vector<int>::const_iterator first = coor_list.begin() + i;
            vector<int>::const_iterator last = coor_list.begin() + i + MAX_IDX;
            vector<int> sub_coor(first, last);
//            cout << " ^^^^^^^^^^^ " << i << ":" << seq_str.size() << endl;
            if(i >= seq_str.size())
                break;
            submakeST(seq_str.substr(i,MAX_IDX), sub_coor);
        }
        for(i = seq_str.size() - MAX_IDX ; i < seq_str.size() - MIN_IDX ; i++){
            vector<int>::const_iterator first = coor_list.begin() + i;
            vector<int>::const_iterator last = coor_list.end();
            vector<int> sub_coor(first, last);
//            cout << " ^^^^^333333^^^^^^ " << i << ":" << seq_str.size() << endl;
            if(i >= seq_str.size())
                break;
            submakeST(seq_str.substr(i), sub_coor);
        }
    }

    bool makeSTtest(SuffixTree &tree, string seq_str, vector<int> coor_list){
            submakeST(seq_str, coor_list);
    }

    void submakeST(string sub_seq, vector<int> sub_coor){
        if(sub_seq.size() < MIN_IDX)
            return;
//        for(int b=0;b<sub_coor.size();b++)
//            cout << sub_coor[b] << " ";
//        cout << sub_seq << endl;

        TreeNode *cur = root;
        int i = 0;
        int j;
        int pivot = -1;
        bool exists_flag;
        map<char, TreeNode*>::iterator it_out;
        map<int, bool>::iterator it;
        while(i < sub_seq.size()) {
//            cout << "8888888" << endl;
//            cout << " 222222 " << endl;
            if(cur->out.size() == 0 || cur->out.count(sub_seq[i]) == 0){
//                cout << " 3332323232333 " << endl;
//                cout << "------ " << sub_seq << " " << sub_seq[i] << " " <<
//                     i << " " << cur->out.count(sub_seq[i]) << endl;
                cur->out[sub_seq[i]] = new TreeNode(sub_seq.substr(i), sub_coor[0]);
                return;
            }
            else if(cur->out.count(sub_seq[i]) == 1 && cur->out[sub_seq[i]]->lab.size() == 1){
//                cout << " 121212211 " << endl;
//#                cur.coor[coor[0]] = True
                cur = cur->out[sub_seq[i]];
                i++;
//                cout << "node " << cur->out[sub_seq[i]]->lab << " : " << cur->out[sub_seq[i]]->lab << endl;
                continue;
            }
            else{
//                cout << " 33333 " << endl;
//                cout << "here " << cur->lab << " : " << cur->out[sub_seq[i]]->lab << endl;
                TreeNode *cur_child = cur->out[sub_seq[i]];

//                for(it_out=cur->out.begin() ; it_out!=cur->out.end() ; ++it_out){
//                    cout << it_out->first << ":" << it_out->second->lab << endl;
//                }
//                cout << sub_seq[i] << ":" << cur->out[sub_seq[i]]->lab << endl;

                exists_flag = true;
                j = 0;
                while(j < cur_child->lab.size() && j < sub_seq.size() - i){
                    if(sub_seq[i+j] != cur_child->lab[j]){
                        exists_flag = false;
                        pivot = j;
//                        cout << "55555555" << " " << i << " " << j << endl ;
                        break;
                    }
                    j = j + 1;
                }
//                cout << "11111 " << cur->out[sub_seq[i]]->lab << " " << cur->out[sub_seq[i]]->lab.size() << " " << j << endl;
                if(exists_flag == true){
//                    cout << " 8888888 " << endl;
//                    cout << "here1111 " << cur->lab << " : " << cur->out[sub_seq[i]]->lab << endl;
//                    if( cur_child->coor.find(sub_coor[0]) == cur_child->coor.end() )
                    cur_child->coor[sub_coor[0]] = true;
                    if(cur_child->lab.size() > j && sub_seq.size() - i > j){
                        cur = cur_child;
                        i += j;
                        continue;
                    }
                    else if(sub_seq.size() - i == j)
                        return;
                    else if(sub_seq.size() - i > j){
                        cur = cur_child;
                        i += j;
                        continue;
                    }
                    else{
                        cur->coor[sub_coor[0]] = true;
                        cout << "should not reach this point" << endl; //, cur_child.lab, seq, i, j)
                        return;
                    }
                }
                else{
                    if(pivot == 0){
                        cout << "wrong : " << sub_seq << " " << cur_child->lab << " " << i << " " << pivot << endl;
//                             << endl << "end_node " << sub_seq.substr(i,i+pivot)
//                                << endl << "start_node " << cur_child->lab.substr(0,pivot) << endl << endl;
                        return;
                    }
//                    cout << " 55555 " << endl;
//                    cout << "here5555 : " << sub_seq << " ::: " << cur_child->lab << " i:" << i << " pivot: " << pivot
//                         << endl << "end_node " << sub_seq.substr(i+pivot)
//                         << endl << "start_node " << cur_child->lab.substr(0,pivot) << endl
//                         << "adjust original "<< cur_child->lab << endl << endl;

                    TreeNode* end_node = new TreeNode(sub_seq.substr(i+pivot), sub_coor[0]); //create end node (new sequence)
                    TreeNode* start_node = new TreeNode(cur_child->lab.substr(0,pivot), sub_coor[0]); //create start node (original sequence)
                    for (it=cur_child->coor.begin() ; it != cur_child->coor.end() ; ++it){
//                        cout << "here22222 " << it->first << endl; // << " : " << cur->out[sub_seq[i]]->lab << endl;
                        if(it->first != 's')
                            start_node->coor[it->first] = true;
                    }
                    start_node->out[sub_seq[i+pivot]] = end_node;
                    start_node->out[cur_child->lab[pivot]] = cur_child;
                    cur_child->lab = cur_child->lab.substr(pivot); //original child's label is curtailed
//#                    cur_child.coor[idx] = True
                    cur->out[sub_seq[i]] = start_node;

//                    cout << "here5555 : " << sub_seq << " ::: " << cur_child->lab << " i:" << i << " pivot: " << pivot
//                         << endl << "end_node " << end_node->lab
//                         << endl << "start_node " << start_node->lab << endl
//                            << "adjust original "<< cur_child->lab << endl << endl;

                    break;
                }

            }

        }
    }

    bool searchPath(string search_seq){
        TreeNode *cur = root;
        int i = 0;
        while(i < search_seq.size()){
            if (cur->out.find(search_seq[i]) == cur->out.end())
                return false; //, []
            else if(cur->out[search_seq[i]]->lab.size() == 1){
                cur = cur->out[search_seq[i]];
                i += 1;
                continue;
            }
            else{
                TreeNode *cur_child = cur->out[search_seq[i]];
                int j = 0;
                while(j < cur_child->lab.size() && j < search_seq.size()-i){
                    if(search_seq[i+j] != cur_child->lab[j])
                        return false; //, []
                    j += 1;
                }
                if(cur_child->lab.size() > j && search_seq.size()-i > j){
                    cur = cur_child;
                    i += j;
                    continue;
                }
                else if(search_seq.size()-i == j){
                    return true; //, cur_child.coor
                }
                else if(search_seq.size()-i > j){
                    cur = cur_child;
                    i += j;
                    continue;
                }
                else{
                    cout << "should not reach this point"; //, cur_child.lab, seq, i, j)
                    return false; //, []
                }
            }
        }
        return true; //, cur_child.coor
    }

//    bool retrievePath(TreeNode &cur_node, string cur_seq){
//        if(cur_node.out.size() == 0){
//            self.printTreeFile.write( cur_seq + cur_node.lab + ':' + ','.join(map(str, cur_node.coor)) + '\n' )
//            return true;
//        }
//        map<char, TreeNode*>::iterator it = cur_node.out.begin();
//        while (it != cur_node.out.end()){
//            if(it->first != 's'){
//                string pass_str = cur_seq;
//                pass_str += cur_node.lab;
//                pass_str += " -> ";
//                retrievePath(*cur_node.out[it->first], pass_str);
//            }
//            else{
//                string pass_str = cur_seq;
//                pass_str += cur_node.lab;
//                pass_str += ":";
//                pass_str += cur_node.coor
//                pass_str += " -> ";
//                retrievePath(*cur_node.out[it->first], pass_str);
//            }
//        }
//    }
//    bool printTree(){
//        TreeNode *cur = root;
//        retrievePath(*cur, "") ;
//    }
};










