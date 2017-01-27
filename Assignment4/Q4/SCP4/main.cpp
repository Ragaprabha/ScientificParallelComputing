#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <iomanip>

using namespace std;

const int m = 4;
const int n = 4;
double **A;

int left(int i, int j);
int right(int i, int j);
int top(int i, int j);
int bottom(int i, int j);
int gaussSeidel(int mynode, int numnodes, int N, double **A, double *x, double *b, double abstol);
double ** CreateMatrix(int rows_local, int N);
void DestroyMatrix(double **A, int rows_local,int N);

int main(int argc, char *argv[]){
    struct timeval start, stop;
    
    int i, j, position, size = (m-1)*(n-1);
    int mynode, totalnodes;
    double w[m+1][n+1];
    double k = 2.0/double(m), h = 1.0/double(n);
    double kSquare = k*k, hSquare = h*h, diagSquare = 2 * kSquare * hSquare;
    double p[size], *b, *x, X[size];
    
    A = new double* [size];
    for (int i=0; i < size; i++)
        A[i] = new double[size];
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    
    gettimeofday(&start, NULL);
    
    for(i = 0; i < size; i++){
        X[i] = 1.0;
    }
    x = X;
    
    for(i = 0; i <= m; i++){
        for(int j=0; j <= n; j++){
            w[i][j] = 0;
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            A[i][j] = 0;
        }
    }
    
    for(j = 1; j < n; j++){
        w[0][j] = exp(j*k);
    }
    
    for(j = 1; j < n; j++){
        w[m][j] = 1.0;
    }
    
    for(i = 1; i < m; i++){
        w[i][0] = 1.0;
    }
    
    for(i = 1; i < m; i++){
        w[i][n] = exp2(2*i*h);
    }
    
    int temp = 0;
    for(j = n-1; j > 0; j--) {
        for(i = 1; i < m; i++){
            double temp1 = ((kSquare * w[i+1][j]) + (kSquare * w[i+1][j]));
            double temp2 = ((hSquare * w[i][j+1]) + (hSquare * w[i][j-1]));
            double temp3 = (2 * (kSquare + hSquare) * w[i][j]);
            double temp4 = (double(i) * k) * (double(i) * k);
            double temp5 = (double(j) * h) * (double(j) * h);
            double temp6 = (kSquare * hSquare * (temp4 + temp5) * exp(double(i) * k * double(j) * h));
            p[temp] = temp1 + temp2 - temp3 + temp6;
            temp++;
        }
    }
    
    b = p;
    position = 0;
    for(j=n-1;j>0;j--){
        for(i=1;i<m;i++){
            int l1=0, l2=0, l3=0, l4=0;
            if((i-1)!= 0){
                l1 =  left(i-1,j);
                A[position][l1-1] = hSquare;
            }
            if((i+1)!= m){
                l2 =  right(i+1,j);
                A[position][l2-1] = hSquare;
            }
            if((j+1)!= n){
                l3 =  top(i,j+1);
                A[position][l3-1] = kSquare;
            }
            if((j-1)!= 0){
                l4 =  bottom(i,j-1);
                A[position][l4-1] = kSquare;
            }
            position++;
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            if(i == j){
                A[i][j] = (-1) * diagSquare;
            }
        }
    }
    int m = gaussSeidel(mynode, totalnodes, size, A, x, b, 1.0);
    gettimeofday(&stop, NULL);
    printf("Time Taken Guass - Seidel (in microseconds):  %lu\n", stop.tv_usec - start.tv_usec);
    MPI_Finalize();

}

int left(int i, int j){
    int l = i+((m-1-j)*(n-1));
    return l;
}

int right(int i, int j){
    int r = i+((m-1-j)*(n-1));
    return r;
}

int top(int i, int j){
    int t = i+((m-1-j)*(n-1));
    return t;
}

int bottom(int i, int j){
    int b = i+((m-1-j)*(n-1));
    return b;
}

int gaussSeidel(int mynode, int numnodes, int N, double **A, double *x, double *b, double abstol){
    int i,j,k,i_global;
    int maxit = 100000;
    int rows_local,local_offset,last_rows_local,*count,*displacements;
    double sum1,sum2,*xold, nomsVector[N], square = 0.0, nomsValue = 0.0;
    double error_sum_local, error_sum_global;
    MPI_Status status;
    
    rows_local = (int) floor((double)N/numnodes);
    local_offset = mynode*rows_local;
    if(mynode == (numnodes-1))
        rows_local = N - rows_local*(numnodes-1);
    
    if(mynode == 0){
        for(i=1;i<numnodes-1;i++){
            for(j=0;j<rows_local;j++)
                MPI_Send(A[i*rows_local+j],N,MPI_DOUBLE,i,j,MPI_COMM_WORLD);
            MPI_Send(b+i*rows_local,rows_local,MPI_DOUBLE,i,rows_local,
                     MPI_COMM_WORLD);
        }
        last_rows_local = N-rows_local*(numnodes-1);
        for(j=0;j<last_rows_local;j++)
            MPI_Send(A[(numnodes-1)*rows_local+j],N,MPI_DOUBLE,numnodes-1,j,
                     MPI_COMM_WORLD);
        MPI_Send(b+(numnodes-1)*rows_local,last_rows_local,MPI_DOUBLE,numnodes-1,
                 last_rows_local,MPI_COMM_WORLD);
    }
    else{
        A = CreateMatrix(rows_local,N);
        x = new double[rows_local];
        b = new double[rows_local];
        for(i=0;i<rows_local;i++)
            MPI_Recv(A[i],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD,&status);
        MPI_Recv(b,rows_local,MPI_DOUBLE,0,rows_local,MPI_COMM_WORLD,&status);
    }
    
    
    xold = new double[N];
    count = new int[numnodes];
    displacements = new int[numnodes];
    
    for(i=0; i<N; i++){
        xold[i] = 1.0;
    }
    
    for(i=0;i<numnodes;i++){
        count[i] = (int) floor((double)N/numnodes);
        displacements[i] = i*count[i];
    }
    count[numnodes-1] = N - ((int)floor((double)N/numnodes))*(numnodes-1);
    
    for(k=0; k<maxit; k++){
        error_sum_local = 0.0;
        for(i = 0; i<rows_local; i++){
            i_global = local_offset+i;
            sum1 = 0.0; sum2 = 0.0;
            for(j=0; j < i_global; j++)
                sum1 = sum1 + A[i][j]*xold[j];
            for(j=i_global+1; j < N; j++)
                sum2 = sum2 + A[i][j]*xold[j];
            
            x[i] = (-sum1 - sum2 + b[i])/A[i][i_global];
            error_sum_local += (x[i]-xold[i_global])*(x[i]-xold[i_global]);
        }
        
        MPI_Allreduce(&error_sum_local,&error_sum_global,1,MPI_DOUBLE,
                      MPI_SUM,MPI_COMM_WORLD);
        MPI_Allgatherv(x,rows_local,MPI_DOUBLE,xold,count,displacements,
                       MPI_DOUBLE,MPI_COMM_WORLD);
        
        for(i=0; i <= (size-1); i++){
            nomsVector[i] = x[i] - xold[i];
        }
        for(int i=0; i <= (size-1); i++){
            square = square + (nomsVector[i] * nomsVector[i]);
        }
        nomsValue = sqrt (square);
        //printf("\nVector Norms: %.2f", nomsValue);
        
        if(sqrt(error_sum_global)<abstol){
            if(mynode == 0){
                for(i=0;i<N;i++)
                    x[i] = xold[i];
            }
            else{
                DestroyMatrix(A,rows_local,N);
                delete[] x;
                delete[] b;
            }
            delete[] xold;
            delete[] count;
            delete[] displacements;
            return k;
        }
    }
    
    if(mynode == 0){
        for(i=0;i<N;i++)
            x[i] = xold[i];
    }
    else{
        DestroyMatrix(A,rows_local,N);
        delete[] x;
        delete[] b;
    }
    delete[] xold;
    delete[] count;
    delete[] displacements;
    
    return maxit;
}

double ** CreateMatrix(int rows_local, int N){
    double **tempA;
    tempA = new double* [rows_local];
    for (int i=0; i < rows_local; i++)
        tempA[i] = new double[N];
    
    for (int i=0; i < rows_local; i++)
    {
        for (int j=0; j < N; j++)
        {
            tempA[i][j] = A[i][j];
        }
    }
    return tempA;
}

void DestroyMatrix(double **A, int rows_local, int N){
    for (int i=0; i < rows_local; i++)
    {
        for (int j=0; j < N; j++)
        {
            A[i][j] = 0.0;
        }
    }
}