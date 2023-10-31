#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define EPS8 1e-8

// Common Function
void printMatrix(double *A, int row, int col);  // Print Matrix
void computeMatrixProduct(double* in1,double* in2, int n1, int n2, int n3, double* out);

void generate_Matrix(int n, double* A); // The matrix in this problem
void SOR_METHOD(int n, double* A, double* b, double* in_x, double omega, double* out_x);

// For Homework3 Problem 6
int main(){
    int n=5;
    double* A=(double*)malloc(sizeof(double)*n*n);
    generate_Matrix(n,A);
    printMatrix(A,n,n);
    free(A);
}




void printMatrix(double *A, int row, int col){
    for(unsigned int i=0;i<row;i++){
        for(unsigned int j=0;j<col;j++){
            if(fabs(A[i*col+j])<EPS8){
                printf("          ");
            }else{
                printf("%7.3f   ",A[i*col+j]);
            }
        }
        printf("\n");
    }
    printf("\n");
}
void computeMatrixProduct(double* in1,double* in2, int n1, int n2, int n3, double* out){
    double sum;
    memset(out, 0, sizeof(double)*n1*n3);
    for(int i=0;i<n1;i++){
        for(int j=0;j<n3;j++){
            for(int k=0;k<n2;k++){
                out[i*n3+j]+=in1[i*n2+k]*in2[k*n3+j];
            }
        }
    }
}
void generate_Matrix(int n, double* A){
    double dn2 = (double)(n*n);
    memset(A, 0, sizeof(double)*n*n);
    for(int i=0;i<n;i++){
        A[i*n+i]=dn2+dn2;
    }
    for(int i=1;i<n-1;i++){
        A[i*n+i+1]=-dn2;
        A[i*n+i-1]=-dn2;
    }
    A[1]=-dn2;
    A[n*n-2]=-dn2;
}

void SOR_METHOD(int n, double* A, double* b, double* in_x, double omega, double* out_x){
    double sum;
    int j;
    for(int i=0;i<n;i++){
        sum=0.0;
        for(j=0;j<i;j++){
            sum += A[i*n+j]*out_x[j];
        }
        for(j=i+1;j<n;j++){
            sum += A[i*n+j]*in_x[j];
        }
        out_x[i] = (1-omega)*in_x[j]+omega*(b[j]-sum)/A[i*n+i];
    }
}
