#include <iostream>
#include <map>
#include <string>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

using namespace std;

struct __align__(8) node{
    char val; 
    struct node *child;
    //map<int, int> comp_map;
	//int sum;
};

__global__
void kernel(node * tree, char *out, int *output_sum, int n)
{
    node *p = tree;
    int i=0;
    while(p->val != 0) {
    out[i++] = p->val;
    //output_sum += p->comp_map[i];
    p = p->child;
    }
}

int main(void)
{
    const int n = 15;
    char data[n] = "tietamattomana";
    node tree[n]; 

	thrust::host_vector<char> h_vec(100);
	//for(int i=0; i<100; i++){
		h_vec[0] = 'R';
		h_vec[1] = 'T';
	//}
	thrust::device_vector<char> d_vec = h_vec;
	thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());

/*
thrust::host_vector<int> h_vec( 16*1024*1024 );
thrust::generate(h_vec.begin(), h_vec.end(), rand);
//  transfer  data  to  the  device
thrust::device_vector<int> d_vec = h_vec;
thrust::sort(d_vec.begin(), d_vec.end());
// sort data on the device
//  transfer  data back  to  host
thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
*/

    node * tree_d;
    char * output_d;
    int * output_sum;
    cudaMalloc((void **)&tree_d, n * sizeof(node));
    cudaMalloc((void **)&output_d, n * sizeof(char));
    cudaMalloc((void **)&output_sum, sizeof(int));

    node * p = tree_d;
    for(int i=0; i<n; i++) {
        tree[i].val = data[i];
        tree[i].child = (++p);
        //tree[i].comp_map[i] = i;
        //tree[i].sum = 0;
    }

    cudaMemcpy(tree_d, tree, n * sizeof(node), cudaMemcpyHostToDevice);
    kernel<<<1,1>>>(tree_d, output_d, output_sum, n);

    char output[n];
    cudaMemcpy(output, output_d, n * sizeof(char), cudaMemcpyDeviceToHost);
    for(int i=0; i<n; i++) {
        std::cout << output[i];
    }
    std::cout << std::endl;

    //int *output_int;
    //cudaMemcpy(output_int, output_sum, sizeof(int), cudaMemcpyDeviceToHost);
    //std::cout << "rrrrrr: " << output_int << std::endl;

    return 0;
}


