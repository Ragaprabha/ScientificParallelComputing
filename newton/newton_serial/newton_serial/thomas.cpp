#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;
double **hessian;
double **hessianInverse;
double *a,*b,*c,*q;
int N;

/*void ThomasAlgorithmLU(int N, double *b, double *a, double *c, double *l, double *u, double *d){
    int i;
    
    
    d[0] = a[0];
    u[0] = c[0];
    for(i=0;i<N-2;i++){
        l[i] = b[i]/d[i];
        d[i+1] = a[i+1] - l[i]*u[i];
        u[i+1] = c[i+1];
    }
    l[N-2] = b[N-2]/d[N-2];
    d[N-1] = a[N-1] - l[N-2]*u[N-2];
    
    return;
}*/

/*double* ThomasAlgorithmSolve(int N, double *l, double *u, double *d, double *x, double *q){
    int i;
    double *y = new double[N];
    
    
    y[0] = q[0];
    for(i=1;i<N;i++)
        y[i] = q[i] - l[i-1]*y[i-1];
    
        x[N-1] = y[N-1]/d[N-1];
    for(i=N-2;i>=0;i--)
        x[i] = (y[i] - u[i]*x[i+1])/d[i];
    
    
    
    delete[] y;
    return x;
}*/



/*double* ThomasAlgorithm(int N, double *b, double *a, double *c, double *x, double *q){
    double *l,*u,*d;
    double *t;
    
    l = new double[N];
    u = new double[N];
    d = new double[N];
    t = new double[N];
    
    ThomasAlgorithmLU(N,b,a,c,l,u,d);
    t = ThomasAlgorithmSolve(N,l,u,d,x,q);
    
    delete[] l;
    delete[] u;
    delete[] d;
    return t;
}*/


int main(int argc, const char * argv[])
{
    printf("Hello, World!\n");
    int n=3;
    N = n+1;
    
    a = new double[N];
    b = new double[N];
    c = new double[N];
    q = new double[N];
    
    
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


    /*printf("Printing the unknown vector before thomas call");
    for(i=0; i<=n;i++)
    {
        printf("%f\n",t[i]);
    }
    
    t = ThomasAlgorithm(N,b,a,c,t,q);
    
    printf("Printing the unknown vector");
    for(i=0; i<=n;i++)
    {
        printf("%f\n",t[i]);
    }*/
        
return 0;

}


