#ifndef _JACOBI_H
#define _JACOBI_H

void jacobi(float **a, int n, float d[], float **v, int *nrot); 
void nrerror(char error_text[]); 
float *vector(long nl, long nh); 
void free_vector(float *v, long nl, long nh);

#endif
