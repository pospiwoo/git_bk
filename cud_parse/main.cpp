#include <iostream>
#include "Graph_modules.h"
#include "Assembly_modules.h"
#include "Indexes.h"
using namespace std;

int main() {
    Indexes *index_obj = Indexes::getInstance();
    index_obj->buildIndexes();
//    index_obj->print_target_index();

    clock_t begin = clock();

    Assembly *assm_obj = new Assembly();
    assm_obj->process();

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time elapsed: " << elapsed_secs << endl;
    return 0;
}









