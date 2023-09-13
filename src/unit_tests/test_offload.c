#include <omp.h>
#include <stdio.h>

#pragma omp declare target
int go_do_something(void);
#pragma omp end declare target

int main ( void ){

  int n;

  printf("Doing something!\n");

#pragma omp target device(0)
  n=go_do_something();


  printf("Did something (%d)!\n",n);
  return 0;
}

int go_do_something(void){

  const int m = 1000;
  int i,j,k[m];

#pragma omp parallel for private(j)
  for (i=0;i<m;i++) {
    j=i*2;
    k[i] = j;
  }

  return k[m-1];
}
