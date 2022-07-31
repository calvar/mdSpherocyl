#include <stdio.h>
#include <stddef.h> 
#include <stdlib.h>
#include <math.h> 

#include "jacobi.h"

#define NR_END 1 
#define FREE_ARG char*
#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau); 


void jacobi(float **a, int n, float d[], float **v, int *nrot) //Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On output,elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a matrix whose columns contain, on output,the normalized eigenvectors of a. nrot returns the number of Jacobi rotations that were required. 
{ 
  int j,iq,ip,i; 
  float tresh,theta,tau,t,sm,s,h,g,c,*b,*z; 
  b=vector(1,n); 
  z=vector(1,n); 
  for (ip=1;ip<=n;ip++) { //Initialize to the ide tity matrix. 
    for (iq=1;iq<=n;iq++) 
      v[ip][iq]=0.0; v[ip][ip]=1.0; 
  } 
  
  for (ip=1;ip<=n;ip++) { //Initialize b and d to the diagonal of a.
    b[ip]=d[ip]=a[ip][ip]; 
    z[ip]=0.0; //This vector will accumulate terms of the form ta pq as i equa- tion (11.1.14). 
  } 
  *nrot=0; 
  for (i=1;i<=50;i++) { 
    sm=0.0; 
    for (ip=1;ip<=n-1;ip++) { //Sum off diagonal elements
      for (iq=ip+1;iq<=n;iq++) 
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) { //The ormal retur ,which relies on quadratic convergence to machi e under  ow. 
      free_vector(z,1,n); 
      free_vector(b,1,n); 
      return; 
    } 
    if (i < 4) 
      tresh=0.2*sm/(n*n); //...on the  rst three sweeps. 
    else 
      tresh=0.0; //...thereafter.
    for (ip=1;ip<=n-1;ip++) { 
      for (iq=ip+1;iq<=n;iq++) { 
	g=100.0*fabs(a[ip][iq]); 
	//After four sweeps,skip the rotatio if the off diagonal element is small 
	if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) 
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) { 
	  h=d[iq]-d[ip]; 
	  if ((float)(fabs(h)+g) == (float)fabs(h)) 
	    t=(a[ip][iq])/h;  //t =1 (2 ) 
	  else { 
	    theta=0.5*h/(a[ip][iq]); //Equation (11.1.10). 
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta)); 
	    if (theta < 0.0) 
	      t = -t; 
	  } 
	  c=1.0/sqrt(1+t*t); 
	  s=t*c; 
	  tau=s/(1.0+c); 
	  h=t*a[ip][iq]; 
	  z[ip] -= h; 
	  z[iq] += h; 
	  d[ip] -= h; 
	  d[iq] += h;
	  a[ip][iq]=0.0; 
	  for (j=1;j<=ip-1;j++) { //Case of rotations 1 dj<p . 
	    ROTATE(a,j,ip,j,iq) 
	  } 
	  for (j=ip+1;j<=iq-1;j++) { //Case of rotations p<j<q . 
	    ROTATE(a,ip,j,j,iq) 
	  } 
	  for (j=iq+1;j<=n;j++) { //Case of rotations q<j dn . 
	    ROTATE(a,ip,j,iq,j) 
	  } 
	  for (j=1;j<=n;j++) { 
	    ROTATE(v,j,ip,j,iq) 
	  } 
	  ++(*nrot); 
	} 
      } 
    } 
    for (ip=1;ip<=n;ip++) { 
      b[ip] += z[ip]; 
      d[ip]=b[ip]; //Update d with the sum of ta pq 
      z[ip]=0.0; //and reinitialize z 
    } 
  } 
  nrerror("Too many iterations in routine jacobi");
}  

void nrerror(char error_text[]) 
/* Numerical Recipes standard error handler */ 
{
  fprintf(stderr,"Numerical Recipes run-time error...\n"); 
  fprintf(stderr,"%s\n",error_text); 
  fprintf(stderr,"...now exiting to system...\n"); 
  exit(1); 
} 

float *vector(long nl, long nh) 
/* allocate a float vector with subscript range v[nl..nh] */ 
{ 
  float *v; 

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float))); 
  if (!v) nrerror("allocation failure in vector()"); 
  return v-nl+NR_END; 
}

void free_vector(float *v, long nl, long nh) 
/* free a float vector allocated with vector() */ 
{ 
  free((FREE_ARG) (v+nl-NR_END)); 
} 
