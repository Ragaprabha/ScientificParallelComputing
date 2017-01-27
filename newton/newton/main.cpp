//
//  Newton-Parallel.cpp
//  Xcode
//
//  Created by Vinod Myll Mylsamy on 5/1/15.
//  Copyright (c) 2015 Vinod Myll Mylsamy. All rights reserved.
//#include <iostream>

#include <stdio.h>
#include <math.h>

using namespace std;
float **hessian;
float **hessianInverse;
double *a,*b,*c,*q,*x;
int N;
// calculate the cofactor of element (row,col)
int GetMinor(float **src, float **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;
    
    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
    
    return 1;
}


// Calculate the determinant recursively.
double CalcDeterminant( float **mat, int order)
{
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];
    
    // the determinant value
    float det = 0;
    
    // allocate the cofactor matrix
    float **minor;
    minor = new float*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new float[order-1];
    
    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        GetMinor( mat, minor, 0, i , order);
        // the recusion is here!
        
        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
        //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }
    
    // release memory
    for(int i=0;i<order-1;i++)
        delete [] minor[i];
    delete [] minor;
    
    return det;
}

// matrix inversioon
// the result is put in Y
float** MatrixInversion(float **A, int order, float **Y)
{
    hessianInverse = new float* [order];
    for (int i=0; i < order; i++)
        hessianInverse[i] = new float[order];
    
    // get the determinant of a
    double det = 1.0/CalcDeterminant(A,order);
    
    // memory allocation
    float *temp = new float[(order-1)*(order-1)];
    float **minor = new float*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));
    
    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            hessianInverse[i][j] = det*CalcDeterminant(minor,order-1);
            if( (i+j)%2 == 1)
                hessianInverse[i][j] = -hessianInverse[i][j];
        }
    }
    /*printf("Printing the inverse matrix \n");
     for(int i=0;i<order;i++)
     {
     for(int j=0; j< order; j++)
     {
     printf("%f\t",Y[i][j]);
     }
     printf("\n");
     }*/
    // release memory
    //delete [] minor[0];
    delete [] temp;
    delete [] minor;
    return hessianInverse;
    
}
void ThomasAlgorithm_P(int mynode, int numnodes,
                       int N, double *b, double *a, double *c, double *x, double *q){
    int i;
    int rows_local,local_offset;
    double S[2][2],T[2][2],s1tmp,s2tmp;
    double *l,*d,*y;
    MPI_Status status;
    
    l = new double[N];
    d = new double[N];
    y = new double[N];
    
    for(i=0;i<N;i++)
        l[i] = d[i] = y[i] = 0.0;
    
    S[0][0] = S[1][1] = 1.0;
    S[1][0] = S[0][1] = 0.0;
    
    rows_local = (int) floor((double)N/numnodes);
    local_offset = mynode*rows_local;
    
    if(mynode==0){
        s1tmp = a[local_offset]*S[0][0];
        S[1][0] = S[0][0];
        S[1][1] = S[0][1];
        S[0][1] = a[local_offset]*S[0][1];
        S[0][0] = s1tmp;
        for(i=1;i<rows_local;i++){
            s1tmp = a[i+local_offset]*S[0][0] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][0];
            s2tmp = a[i+local_offset]*S[0][1] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
            S[0][1] = s2tmp;
        }
    }
    else{
        for(i=0;i<rows_local;i++){
            s1tmp = a[i+local_offset]*S[0][0] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][0];
            s2tmp = a[i+local_offset]*S[0][1] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
            S[0][1] = s2tmp;
        }
    }
    
    for(i=0; i<=log2(numnodes);i++){
        if(mynode+pow(2.0,i) < numnodes)
            MPI_Send(S,4,MPI_DOUBLE,int(mynode+pow(2.0,i)),0,MPI_COMM_WORLD);
        if(mynode-pow(2.0,i)>=0){
            MPI_Recv(T,4,MPI_DOUBLE,int(mynode-pow(2.0,i)),0,MPI_COMM_WORLD,&status);
            s1tmp = S[0][0]*T[0][0] + S[0][1]*T[1][0];
            S[0][1] = S[0][0]*T[0][1] + S[0][1]*T[1][1];
            S[0][0] = s1tmp;
            s1tmp = S[1][0]*T[0][0] + S[1][1]*T[1][0];
            S[1][1] = S[1][0]*T[0][1] + S[1][1]*T[1][1];
            S[1][0] = s1tmp;
        }
    }
    
    d[local_offset+rows_local-1] = (S[0][0] + S[0][1])/(S[1][0] + S[1][1]);
    if(mynode == 0){
        MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(&d[local_offset-1],1,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD,&status);
        if(mynode != numnodes-1)
            MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD);
    }
    
    
    if(mynode == 0){
        l[0] = 0;
        d[0] = a[0];
        for(i=1;i<rows_local-1;i++){
            l[local_offset+i] = b[local_offset+i-1]/d[local_offset+i-1];
            d[local_offset+i] = a[local_offset+i] - l[local_offset+i]*c[local_offset+i-1];
        }
        l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/d[local_offset+rows_local-2];
    }
    else{
        for(i=0;i<rows_local-1;i++){
            l[local_offset+i] = b[local_offset+i-1]/d[local_offset+i-1];
            d[local_offset+i] = a[local_offset+i] - l[local_offset+i]*c[local_offset+i-1];
        }
        l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/d[local_offset+rows_local-2];
    }
    
    /***********************************************************************************/
    
    if(mynode>0)
        d[local_offset-1] = 0;
    
    double * tmp = new double[N];
    for(i=0;i<N;i++)
        tmp[i] = d[i];
    MPI_Allreduce(tmp,d,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<N;i++)
        tmp[i] = l[i];
    MPI_Allreduce(tmp,l,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    delete[] tmp;
    
    if(mynode ==0){    
        /* Forward Substitution [L][y] = [q] */
        y[0] = q[0];
        for(i=1;i<N;i++)
            y[i] = q[i] - l[i]*y[i-1];
        
        /* Backward Substitution [U][x] = [y] */
        x[N-1] = y[N-1]/d[N-1];
        for(i=N-2;i>=0;i--)
            x[i] = (y[i] - c[i]*x[i+1])/d[i];
        
    }
    
    delete[] l;
    delete[] y;
    delete[] d;
    return;
}

int main(int argc, const char * argv[])
{
    printf("Hello, World!\n");
    int n=3;
    N = n+1;
    int totalnodes,mynode;
    

    a = new double[n];
    b = new double[n];
    c = new double[n];
    q = new double[n];
    x = new double[n];
    

    int i,j,k;
    hessian= new float* [n+1];
    for (int i=0; i < n+1; i++)
        hessian[i] = new float[n+1];
    float f[n+1][n+1];
    

    float fn[n+1];
    float x[10]= {1,2,3,4,5,6,7,8,9,10};
    float gradient[n+1];
    
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

    }
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    ThomasAlgorithm_P(mynode,totalnodes,N,b,a,c,x,q);

    if(mynode==0)
    for(int i=0;i<N;i++)
    printf("%f\n",x[i]);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] q;
    delete[] x;
    MPI_Finalize();


    /*hessianInverse = MatrixInversion(hessian,(n+1), hessianInverse);
    printf("Printing the HessianInverse Matrix \n");
    for(i=0; i<=n; i++)
    {
        for(j=0;j<=n;j++)
        {
            printf("%f\t",hessianInverse[i][j]);
        }
        printf("\n");
    }*/
    
    return 0;
    
}


