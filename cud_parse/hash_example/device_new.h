#include <stdio.h>
#define DSIZE 32
__global__ void kernel(){

  int *a = new int[DSIZE];
  for (int i = 0; i < DSIZE; i++) a[i] = i;
  for (int i = 0; i < DSIZE; i++) printf("%d ", a[i]);
  printf("\n");
}

int main(){

  kernel<<<1,1>>>();
  cudaDeviceSynchronize();
}
$ nvcc -o t716 t716.cu
$ cuda-memcheck ./t716
========= CUDA-MEMCHECK
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
========= ERROR SUMMARY: 0 errors
$
