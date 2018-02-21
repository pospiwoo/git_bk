#include <stdio.h>

#define N 2
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

using namespace std;

class vecarray{
    public:
        int *vecptr[N];                //array of pointers pointing to array
        int dim[N];                     //store length of each array pointed to

        __device__ __host__ vecarray(); //constructor
        __device__ __host__ int sum();  //sum up all the elements in the array being
                                       //pointed to
};

vecarray::vecarray(){
    for(int i = 0; i<N; i++)
    {
        vecptr[i] = NULL;
        dim[i] = 0;
    }
}

__device__ __host__ int vecarray::sum(){
    int i=0, j=0, s=0;
    for (i=0; i<N; i++)
        for(j=0; j < dim[i]; j++)
            s += vecptr[i][j];
    return s;
}

__global__ void addvecarray( vecarray * v, int *s){
    *s = v->sum();
}

int main(){                                 //copy *V to device, do sum() and pass back
    vecarray *v, *dev_v;                    //the result by dev_v
    v = new vecarray;
    int a[3] = {1,2,3};                     //initialize v manually
    int b[4] = {4,5,6,7};
    int result = 0;
    int *dev_result;
    v->vecptr[0] = a;
    v->vecptr[1] = b;
    v->dim[0] = 3; v->dim[1] = 4;
    int *vptr[N];

    cudaMalloc((void**)&dev_v, sizeof(vecarray));
    cudaCheckErrors("cudaMalloc1 fail");
    cudaMemcpy(dev_v, v, sizeof(vecarray),cudaMemcpyHostToDevice); //copy class object
    cudaCheckErrors("cudaMemcpy1 fail");

    for(int i = 0; i < N; i++){
        cudaMalloc((void**)&(vptr[i]), v->dim[i]*sizeof(int));
        cudaCheckErrors("cudaMalloc2 fail");
        cudaMemcpy(&(dev_v->vecptr[i]), &vptr[i], sizeof(int*), cudaMemcpyHostToDevice);
        cudaCheckErrors("cudaMemcpy2 fail");
    }

    for(int i = 0; i<N; i++ ){                   //copy arrays
        cudaMemcpy(vptr[i], v->vecptr[i], v->dim[i]*sizeof(int), cudaMemcpyHostToDevice);
        cudaCheckErrors("cudaMemcpy3 fail");
    }
    cudaMalloc((void **)&dev_result, sizeof(int));
    cudaCheckErrors("cudaMalloc3 fail");
    addvecarray<<<1,1>>>(dev_v, dev_result);

    cudaMemcpy(&result, dev_result, sizeof(int), cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy4 fail");
    printf("the result is %d\n", result);
    return 0;
}

