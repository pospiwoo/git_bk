#include <string>
#include "Gene.h"
using namespace std;

int main() {
    ConstructGraph *ConstructGraph_inst = new ConstructGraph;
    cout << "Parsing DNA" << endl;
    ConstructGraph_inst->parseDNA();
    //ConstructGraph_inst->referenceG();
    //ConstructGraph_inst->testReferenceG();
    cout << "Build GGraph DNA" << endl;
    ConstructGraph_inst->parseGFF();
//    ConstructGraph_inst->parseGFFtest();
    return 0;
}
