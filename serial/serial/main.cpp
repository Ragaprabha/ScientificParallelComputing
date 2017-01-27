//
//  main.cpp
//  serial
//
//  Created by Vinod Myll Mylsamy on 5/6/15.
//  Copyright (c) 2015 Vinod Myll Mylsamy. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Utilities.cuh"

using namespace std;
double **hessian;
double **hessianInverse;
double *a,*b,*c,*q,*t;
double *d_a,*d_b,*d_c,*d_l,*d_u,*d_d;
int N;

// Device code
_global_ void ThomasAlgorithmLU(int N, double *b, double *a, double *c, double *l, double *u, double *d){
    
    int i = blockDim.x * blockIdx.x + threadIdx.x;
 
 
    d[0] = a[0];
    u[0] = c[0];
    if(i< N-2){
     l[i] = b[i]/d[i];
     d[i+1] = a[i+1] - l[i]*u[i];
     u[i+1] = c[i+1];
    }
    
    l[N-2] = b[N-2]/d[N-2];
    d[N-1] = a[N-1] - l[N-2]*u[N-2];
 
    return;
 }

void ThomasAlgorithmSolve(int N, double *l, double *u, double *d, double *t, double *q){
 int i;
 double *y = new double[N];
 
 
 y[0] = q[0];
 for(i=1;i<N;i++)
 y[i] = q[i] - l[i-1]*y[i-1];
 
 t[N-1] = y[N-1]/d[N-1];
 for(i=N-2;i>=0;i--)
 t[i] = (y[i] - u[i]*t[i+1])/d[i];
 
 
 
 delete[] y;
 return;
 }



void ThomasAlgorithm(int N, double *b, double *a, double *c, double *t, double *q){
 double *l,*u,*d;
    int i;
 int size =(int) N * sizeof(double);
 
 l = new double[N];
 u = new double[N];
 d = new double[N];
    
    for(i=0;i<N;i++)
    {
        l[i] = 0;
    }
    for(i=0;i<N;i++)
    {
        u[i] = 0;
    }
    for(i=0;i<N;i++)
    {
        d[i] = 0;
    }
    
     // Alloc space for device copies
    cudaMalloc((void **)&d_a, size);
    cudaMalloc((void **)&d_b, size);
    cudaMalloc((void **)&d_c, size);
    cudaMalloc((void **)&d_l, size);
    cudaMalloc((void **)&d_u, size);
    cudaMalloc((void **)&d_d, size);
    
    // Copy inputs to device
    cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_c, c, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_l, q, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, t, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_d, N, size, cudaMemcpyHostToDevice);
    
    //Invoke kernel
    int threadsPerBlock = 2;
    int blocksPerGrid = (N + threadsPerBlock -1)/ threadsPerBlock;
    
    //Invoke Kernel
    ThomasAlgorithmLU<<<blocksPerGrid,threadsPerBlock>>>(N,b,a,c,l,u,d);
    
    
    // Copy result back to host
    cudaMemcpy(l, d_l, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(d, d_d, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost);
    
    
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    cudaFree(d_l);
    cudaFree(d_u);
    cudaFree(d_d);
    
    ThomasAlgorithmSolve(N,l,u,d,t,q);
 
    delete[] l;
    delete[] u;
    delete[] d;
    
    return;
 }


int main(int argc, const char * argv[])
{
    printf("Hello, World!\n");
    int n=3;
    N = n+1;
    
    a = new double[N];
    b = new double[N];
    c = new double[N];
    q = new double[N];
    t = new double[N];
    
    int i,j;
    hessian= new double*[n+1];
    for (int i=0; i < n+1; i++)
        hessian[i] = new double[n+1];
    double f[n+1][n+1];
    
    
    double fn[n+1];
    double x[10]= {1,2,3,4,5,6,7,8,9,10};
    double gradient[n+1];
    
    for (i=0;i<n;i++)
    {
        if(i==0)
        {
            fn[i]= (100 * ((x[i]*x[i])- x[i+1]) * ((x[i]*x[i])- x[i+1]))+ ((x[i]-1) * (x[i]-1));
        }
        else
        {
            fn[i] = fn[i-1] + (100 * ((x[i]*x[i])- x[i+1]) * ((x[i]*x[i])- x[i+1]))+ ((x[i]-1) * (x[i]-1));
        }
    }
    for(i=0;i<=n;i++)
    {
        for(j=0;j<=n;j++)
        {
            if(i==0 && j==0)
            {
                f[i][j] = (800 * x[0]*x[0])+ (400 * ((x[0]*x[0])-x[1]))+2;
            }
            else if(i== n && j== n)
            {
                f[i][j] = 200;
            }
            else if (i==j && i!=n)
            {
                f[i][j] = (400 *((x[i]*x[i])-x[i+1]))+ (800* x[i]*x[i])+202;
            }
            else if (j== (i-1)|| j == (i+1))
            {
                f[i][j] = -400 * x[i];
            }
            else{
                f[i][j] = 0;
            }
        }
    }
    for(i=0;i<=n;i++)
    {
        if(i==0)
        {
            gradient[0] = (400 * x[0] * ((x[0] * x[0])-x[1])) + (2* (x[0] -1 ));
        }
        else if (i== n)
        {
            gradient[n] = (-200 * (( x[n-1] * x[n-1])- x[n]));
        }
        else if ( i < n)
        {
            gradient[i] = (-200 * ((x[i-1] * x[i-1]) - x[i])) + ( 400 * x[i] * (( x[i] * x[i])-x[i+1])) + ( 2* (x[i] - 1));
        }
    }
    printf("Printing the gradient vector \n");
    for (i=0; i<=n;i++)
    {
        printf("%f\n",gradient[i]);
    }
    printf("Printing the HessianMatrix\n");
    for(i=0; i<=n; i++)
    {
        for(j=0;j<=n;j++)
        {
            hessian[i][j] = f[i][j];
            printf("%f\t",hessian[i][j]);
        }
        printf("\n");
    }
    for(i=0; i <= n; i++)
    {
        j=i;
        a[i]= hessian[i][j];
        
    }
    printf("Printing the diagonal vector");
    for(i=0;i<=n;i++)
    {
        printf("%f\n",a[i]);
    }
    b[0]=0;
    for(i=1; i <= n; i++)
    {   j= i-1;
        b[i]= hessian[i][j];
        
    }
    printf("Printing the lower diagonal vector");
    for(i=0;i<=n;i++)
    {
        printf("%f\n",b[i]);
    }
    for(i=0; i < n; i++)
    {   j= i+1;
        c[i]= hessian[i][j];
        
    }
    c[n]=0;
    printf("Printing the upper diagonal vector");
    for(i=0;i<=n;i++)
    {
        printf("%f\n",c[i]);
    }
    for(i=0;i<=n;i++)
    {
        q[i] = - gradient[i];
    }
    printf("Printing the right hand side vector");
    for(i=0;i<=n;i++)
    {
        printf("%f\n",q[i]);
    }
    for(i=0;i<=n;i++)
    {
        t[i] = 0;
    }
    
     printf("Printing the unknown vector before thomas call");
     for(i=0; i<=n;i++)
     {
     printf("%f\n",t[i]);
     }
    
    ThomasAlgorithm(N,b,a,c,t,q);
    
     printf("Printing the unknown vector");
     for(i=0; i<=n;i++)
     {
     printf("%f\n",t[i]);
     }
    
    // Cleanup
    free(a);
    free(b);
    free(c);
    free(q);
    free(t);
    

    return 0;
    
}


