#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

using namespace std;

int main(int argc, const char * argv[]) {
   
    int i,j,k;
    double x = -1.2,y = 1;
    
    double xnew[2];
    double xold[2]= {-1.2,1};
    double xdiff[2];
    double gradient[2];
    double hessian[2][2];
    double hessianInverse[2][2];
    double hessianInv[2][2];
    double s[2];
    double f;
    
    for(k=0;k<20;k++)
    {
        f = 0;
        f = ((1- x)*(1- x)) + (100*(y-(x*x))*(y-(x*x)));
        printf("The value of the function is %f\n\n",f);
        
        double g1 = (double)((-2.0*(1.0-x))-((400.0*x)*(y-(x*x))));
        double g2 =(double)(200*(y-(x*x)));
        double h1 = (double)(2-(400*y)+(1200*x*x));
        double h2 = (double)(-400 *x);
        double h3 = (double)(-400 * x);
        double h4 = (double)(200);
        double determinant = 1/((h1 * h4)-(h2*h3));
        
        
        gradient[0] = g1;
        gradient[1] = g2;
        hessian[0][0]=h1;
        hessian[0][1]=h2;
        hessian[1][0]= h3;
        hessian[1][1]=h4;
        hessianInv[0][0]=h4;
        hessianInv[0][1]=-h3;
        hessianInv[1][0]=-h2;
        hessianInv[1][1]=h1;
        printf("Gradient\n");
        for (i=0; i<=1;i++)
        {
            printf("%f\n",gradient[i]);
        }
        printf("\nHessianMatrix\n");
        for(i=0; i<=1; i++)
        {
            for(j=0;j<=1;j++)
            {
                printf("%f\t",hessian[i][j]);
            }
            printf("\n");
        }
        
        for(i=0; i<=1; i++)
        {
            for(j=0;j<=1;j++)
            {
                hessianInverse[i][j] = determinant * hessianInv[i][j];
            }
        }
        for(i=0; i<=1; i++)
        {   s[i] = 0 ;
            for(j=0;j<=1;j++)
            {
                s[i] = s[i]+ (hessianInverse[i][j] * (- gradient[j]));
            }
        }
        printf("\nNewton Step: \n");
        
        for(i=0; i<=1; i++)
        {
            printf("%f\n", s[i]);
        }
        
        printf("\nThe Value of X after %d iteration \n", k);
        for(i=0; i<=1; i++)
        {
            xnew[i] = 0;
            xnew [i] = s[i] + xold [i];
        }
        for(i=0; i<=1; i++)
        {
            printf("%f\n", xnew[i]);
        }
        for(i=0; i<=1; i++)
        {
            xdiff[i] = xnew[i] - xold[i];
        }
        for(i=0; i<=1; i++)
        {
            xold[i] = xnew[i];
        }
        x = xnew[0];
        y = xnew[1];
        
        double square1 = 0, nomsValue1 = 0;
        for(int i=0; i <= 1; i++){
            square1 = square1 + (xdiff[i] * xdiff[i]);
        }
        nomsValue1 = sqrt (square1);
        printf("\nVector Noms: %.2f\n", nomsValue1);
        
        if(nomsValue1 <= 0){
            break;
        }
        
        
    }
    return 0;
}

