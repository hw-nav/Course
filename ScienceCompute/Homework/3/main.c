#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define EPS8 1e-8

// Common Function
void printMatrix(double *A, int row, int col);  // Print Matrix
void copyMatrix(int length, double *A, double *B);
void computeMatrixProduct(int n1, int n2, int n3, double* in1,double* in2, double* out);
void computeMatrixAdd(int length, double* in1, double* in2, double* out);
void computeMatrixMinus(int length, double* in1, double* in2, double* out);
double computeMatrixInfiteNorm(int row,int col, double* A);
int GaussEliminationInverse(int n, double* A,double* invA);

void generate_Matrix(int n, double* A); // The matrix in this problem
void get_Iter_Matrix(int n, double* A, double* B_J, double* B_GS);
void SOR_METHOD(int n, double* A, double* b, double* in_x, double omega, double* out_x);

// For Homework3 Problem 6

// ERROR: g++ -lm main.c
int main(){
    int n;
    double *A, *B_J, *B_GS, *I1, *I2;
    double rho;

    // (c)
    n = 4;
    A=(double*)malloc(sizeof(double)*n*n);
    B_J=(double*)malloc(sizeof(double)*n*n);
    B_GS=(double*)malloc(sizeof(double)*n*n);
    I1=(double*)malloc(sizeof(double)*n*n);
    I2=(double*)malloc(sizeof(double)*n*n);

    generate_Matrix(n,A);
    get_Iter_Matrix(n,A,B_J,B_GS);

    printf("Jacobi:\n");
    copyMatrix(n*n,B_J,I2);
    rho = computeMatrixInfiteNorm(n,n,I2);
    printf(">>> k= 1, rho=%7.3f\n",rho);
    copyMatrix(n*n,I2,I1);
    for(int k=2;k<=20;k++){
        computeMatrixProduct(n,n,n,I1,B_J,I2);
        rho = computeMatrixInfiteNorm(n,n,I2);
        rho = pow(rho,1.0/k);
        printf(">>> k=%2d, rho=%7.3f\n",k,rho);
        copyMatrix(n*n,I2,I1);
    }

    printf("Gauss-Siedel:\n");
    copyMatrix(n*n,B_GS,I2);
    rho = computeMatrixInfiteNorm(n,n,I2);
    printf(">>> k= 1, rho=%7.3f\n",rho);
    copyMatrix(n*n,I2,I1);
    for(int k=2;k<=20;k++){
        computeMatrixProduct(n,n,n,I1,B_GS,I2);
        rho = computeMatrixInfiteNorm(n,n,I2);
        rho = pow(rho,1.0/k);
        printf(">>> k=%2d, rho=%7.3f\n",k,rho);
        copyMatrix(n*n,I2,I1);
    }

    free(A);
    free(B_J);
    free(B_GS);
    free(I1);
    free(I2);
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
void copyMatrix(int length, double *A, double *B){
    for(int i=0;i<length;i++){
        B[i]=A[i];
    }
}
void computeMatrixProduct(int n1, int n2, int n3, double* in1,double* in2, double* out){
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
void computeMatrixAdd(int length, double* in1, double* in2, double* out){
    for(int i=0;i<length;i++){
        out[i]=in1[i]+in2[i];
    }
}
void computeMatrixMinus(int length, double* in1, double* in2, double* out){
    for(int i=0;i<length;i++){
        out[i]=in1[i]-in2[i];
    }
}
double computeMatrixInfiteNorm(int row,int col, double* A){
    double sum,res=0.0;
    for(int i=0;i<row;i++){
        sum=0.0;
        for(int j=0;j<col;j++){
            sum+=A[i*col+j];
        }
        if (res<fabs(sum)){
            res=fabs(sum);
        }
    }
    return res;
}
int GaussEliminationInverse(int n, double* A, double* invA) {
    double* copyA=(double*)malloc(sizeof(double)*n*n); // copy of A
    double div = 0.0, factor;

    // copy
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            copyA[i*n+j]=A[i*n+j];
        }
    }
	
    // Initial invA
	memset(invA, 0, sizeof(double)*n*n);
	for(int i = 0; i < n; i++) {
		invA[i*n + i] = 1.0;
	}

    // A L-trio
	for (int i = 0; i < n; i++) {
        // Row i
        if(fabs(copyA[i*n+i])<1e-15){
            return 1;                   // < 1e-15 error return 1
        }
		div = 1.0 / copyA[i*n + i];

        // Time limit to prioritization
        for(int j=0;j<n;j++){
            copyA[i*n+j] *= div;
            invA[i*n+j] *= div;
        }
		for(int k=i+1;k<n;k++){
            factor = copyA[k*n+i];
            for(int j=0;j<n;j++){
                copyA[k*n+j] -= factor*copyA[i*n+j];
                invA[k*n+j] -= factor*invA[i*n+j];
            }
        }
	}

    // Inverse
    for(int i=n-1;i>=0;i--){
        for(int k=0;k<=i-1;k++){
            factor = copyA[k*n+i];
            for(int j=0;j<n;j++){
                copyA[k*n+j] -= factor*copyA[i*n+j];
                invA[k*n+j]-=factor*invA[i*n+j];
            }
        }
    }

    free(copyA);
	return 0;
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

void get_Iter_Matrix(int n, double* A, double* B_J, double* B_GS){
    double* L = (double*)malloc(sizeof(double)*n*n);
    double* D = (double*)malloc(sizeof(double)*n*n);
    double* U = (double*)malloc(sizeof(double)*n*n);
    double* Inv = (double*)malloc(sizeof(double)*n*n);
    double* L_U = (double*)malloc(sizeof(double)*n*n);
    memset(L, 0, sizeof(double)*n*n);
    memset(D, 0, sizeof(double)*n*n);
    memset(U, 0, sizeof(double)*n*n);
    memset(L_U, 0, sizeof(double)*n*n);

    int i,j;
    for(i=0;i<n;i++){
        D[i*n+i]=A[i*n+i];
        for(j=0;j<i;j++){
            L[i*n+j]=-A[i*n+j];
            U[j*n+i]=-A[j*n+i];
        }
    }
    

    computeMatrixAdd(n*n,L,U,L_U);
    GaussEliminationInverse(n,D,Inv);
    computeMatrixProduct(n,n,n,Inv,L_U,B_J);
    free(L_U);


    double* D_L = (double*)malloc(sizeof(double)*n*n);
    memset(D_L, 0, sizeof(double)*n*n);
    computeMatrixMinus(n*n,D,L,D_L);
    free(L);
    free(D);

    GaussEliminationInverse(n,D_L,Inv);
    free(D_L);
    computeMatrixProduct(n,n,n,Inv,U,B_GS);
    free(U);
    free(Inv);
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