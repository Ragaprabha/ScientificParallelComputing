#include <stdio.h>
#include <iostream>
#include <math.h>
#include <sys/time.h>

using namespace std;

#define PI 3.14159265

const int node = 4;

int main(){
  struct timeval start, stop;

  int i,j,k;
  int index1,index2,offset, block;
  int rStart=0, rEnd=3, cStart = 0, cEnd = 3;
  int size = node * node;
  double h = 1/(double) node;
  double U = h*h, alpha, beta;
  double matrix[size][size], matrix1[size*2][size*2], d;
  double qVec[size], temp[size], r0[size], x0[size], tempx[size], tempr[size], p[size];
  double inverse[size][size], z0[size], tempz[size];

  gettimeofday(&start, NULL);

  /* Matrix Generation */

  for(i=0; i<=(size-1); i++){
    for(j=0; j<=(size-1); j++){
     matrix[i][j] = 0; 
     inverse[i][j] = 0; 
    }
  }
  
  if(node == 4){
  for(block=1; block<=size;block++){
    if(block == 1 || block == 6 || block == 11 || block == 16){ 
      for(i= rStart; i<=rEnd; i++){
        for(j= cStart; j<=cEnd; j++){
          if(i == j){
            matrix[i][j] = (double) node + U;
            matrix[i][j+1] = -1.0;
            if(j != 0){
              matrix[i][j-1] = -1.0;
            }       
          }else if((j == cEnd && i == rStart) || (i == rEnd && j == cStart)){
             matrix[i][j] = -1.0;
          }
        }
      }
    } else if(block == 2 || block == 4 || block == 5 || block == 7 || block == 10 || block == 12 || block == 13 || block == 15){
       for(i= rStart; i<=rEnd; i++){
        for(j=cStart; j<=cEnd; j++){
          if(j == (i+4) || j == (i+12) || i == (j+4) || i == (j+12)){
            matrix[i][j] = -1.0;
          }
        }
      }
    }

    if(cEnd == size-1){
      cStart = 0;
      cEnd = node-1;
      rStart = rStart + node;
      rEnd = rEnd + node;
    }else{
      cStart = cStart + node;
      cEnd = cEnd + node;
    }
  }}else if(node == 20){
    for(i=0; i<(size-1); i++){
      k = 1;
      for(j=0; j<(size-1); j++){
         if(i == j){
            matrix[i][j] = (double) node + U;
            matrix[i][j+1] = -1.0;
            if(j != 0){
              matrix[i][j-1] = -1.0;
            }
         }else if(j == (i+(20*k)) || i == (j+(20*k))){
            matrix[i][j] = -1.0;
            k = k +2;
         }          
      }
    }
  }
  
  for(i=0; i<=(size-1); i++){
    for(j=0; j<=(size-1); j++){
      printf("%.2f\t", matrix[i][j]);   
    }
    printf("\n");
  }

  /*Inverse of a Matrix */

  for(i=0; i < size; i++){
    for(j=0; j< size; j++){
     matrix1[i][j] = matrix[i][j];  
    }
  }
  
  for (i = 0; i < size; i++){
    for (j = size; j < size * 2; j++){
       matrix1[i][j] = 0;
    }
  }
  
  for (i = 0; i < size; i++){
    for (j = 0; j < 2 * size; j++){
      if (j == (i + size)){
        matrix1[i][j] = 1.0;
      }
    }
  }
  
  for (i = size-1; i > 0; i--){
    if (matrix1[i-1][0] < matrix1[i][0]){
      for(j = 0; j < size * 2; j++){
        d = matrix1[i][j];
        matrix1[i][j] = matrix1[i-1][j];
        matrix1[i-1][j] = d;
      }
    }
  }
 
  for (i = 0; i < size; i++){
    for (j = 0; j < size * 2; j++){
      if (j != i){
        d = matrix1[j][i] / matrix1[i][i];
        for (k = 0; k < size * 2; k++){
          matrix1[j][k] = matrix1[j][k] - (matrix1[i][k] * d);
        }
      }
    }
  }
    
  for (i = 0; i < size; i++){
    d=matrix1[i][i];
    for (j = 0; j < size * 2; j++){
      matrix1[i][j] = matrix1[i][j] / d;
    }
  }
  
  for (i = 0; i < size; i++){
    for (j = size; j < size * 2; j++){   
      inverse[i][j-(size-1)]= matrix1[i][j];
    }
  }

  /* Q[i][j] Matrix */
  
  int m = 0, n = 0;
  double tempm = 0, tempn = 0;
  for(i=0; i <= (size-1); i++){
    tempm = sin (2*PI*h*m);
    tempn = sin (2*PI*h*n);
    qVec[i] = ((8*(PI*PI)+1)*(h*h)*tempm*tempn);
    if(m == (node-1)){
      m = 0;
      n = n+1;
    }else{
      m=m+1;
    }  
  }
  
  /* calculate r0 Matrix*/

  double temp1, temp2, temp3, temp4;
  int a;
  for(a=1; a<=size; a++){
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;

    if(a == 1){
      for(i=0; i <= (size-1); i++){
        x0[i] = i;
      }

      for(i=0; i <= (size-1); i++){
        for(j=0; j <= (size-1); j++){
          temp[i] = temp[i] + (matrix[i][j] * x0[j]);
        }
      }
      for(i=0; i <= (size-1); i++){
        r0[i] = qVec[i] - temp[i];
      }
      
      /* Calculate r0 */
     
      for(i=0; i <= (size-1); i++){
        r0[i] = qVec[i] - temp[i];
      }
      
      /* Calculate z0 */
     
      for(i=0; i <= (size-1); i++){
        for(j=0; j <= (size-1); j++){
          z0[i] = z0[i] + (inverse[i][j] * r0[j]);
        }
      }

      for(i=0; i <= (size-1); i++){
        p[i] = z0[i];
      }

      /* Calculate scalar Alpha */
      
      for(i=0; i <= (size-1); i++){
        for(j=0; j <= (size-1); j++){
          temp[i] = temp[i] + (matrix[i][j] * p[j]);
        }
      }
      for(i=0; i <= (size-1); i++){
        temp1 = temp1 + (r0[i] * z0[i]);
        temp2 = temp2 + (p[i] * temp[i]);
      } 
      alpha = temp1/temp2;
      
      for(i=0; i <= (size-1); i++){
         temp[i] = alpha * p[i];
         tempx[i] = x0[i];
         x0[i] = temp[i] + tempx[i];
       }
    }else{    
       for(i=0; i <= (size-1); i++){
         for(j=0; j <= (size-1); j++){
           temp[i] = temp[i] + (matrix[i][j] * p[j]);
         }
       }
       for(i=0; i <= (size-1); i++){
	 temp[i] = alpha * temp[i];
         tempr[i] = r0[i];
         r0[i] = tempr[i] - temp[i];
         tempz[i] = z0[i];
       }

       for(i=0; i <= (size-1); i++){
         for(j=0; j <= (size-1); j++){
            z0[i] = z0[i] + (inverse[i][j] * r0[j]);
        }
      }      

       for(i=0; i <= (size-1); i++){
         temp3 = temp3 + (z0[i] * r0[i]);
         temp4 = temp4 + (tempz[i] * tempr[i]);
       }
       beta = temp3/temp4;
       for(i=0; i <= (size-1); i++){
         temp[i] = beta * p[i];
         p[i] = z0[i] + temp[i];
       }
       for(i=0; i <= (size-1); i++){
         for(j=0; j <= (size-1); j++){
           temp[i] = temp[i] + (matrix[i][j] * p[j]);
         }
       }

       for(i=0; i <= (size-1); i++){
        temp1 = temp1 + (z0[i] * r0[i]);
        temp2 = temp2 + (p[i] * temp[i]);
       } 
       alpha = temp1 / temp2;
       for(i=0; i <= (size-1); i++){
         temp[i] = alpha * p[i];
         tempx[i] = x0[i];
         x0[i] = temp[i] + tempx[i];
       }  
     }    
   } 

  printf("\n");
  for(i=0; i <= (size-1); i++){
    printf("%.2f\n", x0[i]);
  }
  printf("\n");
  gettimeofday(&stop, NULL);
  printf("Time Taken (in microseconds):  %lu\n", stop.tv_usec - start.tv_usec);
}
  