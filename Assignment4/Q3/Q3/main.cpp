#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <iomanip>

using namespace std;const int m = 4;
const int n = 4;
double **A;
double **B;

int left(int i, int j);
int right(int i, int j);
int top(int i, int j);
int bottom(int i, int j);
void ThomasAlgorithm(int N, double b, double a, double c, double *x, double *q);
int main(int argc, char *argv[]){
    struct timeval start, stop;
    
    int i, j, position, size = (m-1)*(n-1);
    double w[m+1][n+1];
    double k = 1.0/double(m), h = 0.25;
    double alpha = 1;
    double kSquare = k*k, hSquare = h*h, diagSquare = 2 * kSquare * hSquare;
    double p[size],r[size], *b, *x, X[size];
    
    A = new double* [size];
    for (int i=0; i < size; i++)
        A[i] = new double[size];
    
    B = new double* [size];
    for (int i=0; i < size; i++)
        B[i] = new double[size];
    
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
            A[i][j] = 0.0;
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            B[i][j] = 0.0;
        }
    }
    
    for(j = 1; j < n; j++){
        w[0][j] = 0;
    }
    
    for(j = 1; j < n; j++){
        w[m][j] = 0;
    }
    
    for(i = 1; i < m; i++){
        w[i][0] = 0;
    }
    
    for(i = 1; i < m; i++){
        w[i][n] = 0;
    }
    
    int temp = 0;
    for(j = n-1; j > 0; j--) {
        for(i = 1; i < m; i++){
            double temp1 = (5*w[i][j+1])+(3*w[i][j]);
            double temp2 = 2*(w[i+1][j]+w[i-1][j]+w[i+1][j+1]+w[i-1][j+1]);
            double temp3 = 0.0;
            if((i*k)<= 0.5){
                temp3 = 1;
            }
            if((i*k)>0.5){
                temp3 =0;
            }
            p[temp] = temp1 - temp2 - (0.5 * temp3);
            temp++;
        }
    }
    
    temp = 0;
    for(j = n-1; j > 0; j--) {
        for(i = 1; i < m; i++){
            double temp1 = (5*w[i][j+1])+(3*w[i][j]);
            double temp2 = 2*(w[i+1][j]+w[i-1][j]+w[i+1][j+1]+w[i-1][j+1]);
            double temp3 = 0.0;
            if((i*k)<= 0.5){
                temp3 = 0;
            }
            if((i*k)>0.5){
                temp3 =1;
            }
            r   [temp] = temp1 - temp2 - (0.2 * temp3);
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
                A[position][l1-1] = -2.0;
            }
            if((i+1)!= m){
                l2 =  right(i+1,j);
                A[position][l2-1] = -2.0;
            }
            if((j+1)!= n){
                l3 =  top(i,j+1);
                A[position][l3-1] = -2.0;
            }
            if((j-1)!= 0){
                l4 =  bottom(i,j-1);
                A[position][l4-1] = -2.0;
            }
            position++;
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            if(i == j){
                A[i][j] = 3.0;
                A[i][j+1]=5.0;
            }
        }
    }
    
    
    position = 0;
    for(j=n-1;j>0;j--){
        for(i=1;i<m;i++){
            int l1=0, l2=0, l3=0, l4=0;
            if((i-1)!= 0){
                l1 =  left(i-1,j);
                B[position][l1-1] = -2.0;
            }
            if((i+1)!= m){
                l2 =  right(i+1,j);
                B[position][l2-1] = -2.0;
            }
            if((j+1)!= n){
                l3 =  top(i,j+1);
                B[position][l3-1] = -2.0;
            }
            if((j-1)!= 0){
                l4 =  bottom(i,j-1);
                B[position][l4-1] = -2.0;
            }
            position++;
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            if(i == j){
                B[i][j] = 3.0;
                B[i][j+1]=5.0;
            }
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
                A[i][j] = A[i][j]+B[i][j];
        }
    }
    
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            printf("%f\t", A[i][j]);
        }
        printf("\n");
    }
    
    for(i = 0; i < size; i++){
        p[i]=p[i]+r[i];
    }
    printf("\n");
    
    for(i = 0; i < size; i++){
        printf("%f\n", r[i]);
    }
    
    ThomasAlgorithm(size, 6.0, 10,-4.0, x,b);
    
    gettimeofday(&stop, NULL);
    printf("Time Taken Thomas Algorithm (in microseconds):  %lu\n", stop.tv_usec - start.tv_usec);
    
    
}
void ThomasAlgorithm(int N, double b, double a, double c, double *x, double *q){
    int i;
    double *l,*u,*d,*y;
    
    l = new double[N];
    u = new double[N];
    d = new double[N];
    y = new double[N];
    
    u[0] = c;
    for(i=0;i<N-2;i++){
        l[i] = b/d[i];
        d[i+1] = a - l[i]*u[i];
        u[i+1] = c;
    }
    l[N-2] = b/d[N-2];
    d[N-1] = a - l[N-2]*u[N-2];
    
    
    y[0] = q[0];
    for(i=1;i<N;i++)
        y[i] = q[i] - l[i-1]*y[i-1];
    
   
    x[N-1] = y[N-1]/d[N-1];
    for(i=N-2;i>=0;i--)
        x[i] = (y[i] - u[i]*x[i+1])/d[i];
    
    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    return;
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



