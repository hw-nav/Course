#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define EPS8 1e-8
#define MAXITER 100000

void generateA(int n, double* A);
double fun(double x);
void generatef(int n, double* f);
void steepestDescentIter(int n, double* A, double* r, double* x);
void conjugateGradientIter(int n,double* A, double* r, double* rTr_iter, double* d, double* x);
double inftyNorm(int n, double *x);


// ERROR: g++ -lm main.c
int main(){
    //int len=1, nlist[]={5};
    int len=6, nlist[]={9,19,39,79,159,319};
    int n;
    int iter;
    double *A, *f, *u;
    double *r,*d,rTr_iter;
    double norm, norm0;

    // Test with Question 2, iter=153
    // double B[]={2.0,-1.0,0.0,-1.0,3.0,-2.0,0.0,-2.0,2.0};
    // double b[]={3.0,3.0,-4.0};
    // r = (double*)malloc(sizeof(double)*3);
    // d = (double*)malloc(sizeof(double)*3);
    // double x[]={1.0,1.0,1.0};

    // rTr_iter = 0.0;
    // for(int i=0;i<3;i++){
    //     r[i] = b[i] - B[i*3+0]*x[0] - B[i*3+1] * x[1] - B[i*3+2];
    //     rTr_iter += r[i]*r[i];
    //     d[i] = r[i];
    // }
    // norm0 = inftyNorm(3,r);
    // iter=1;
    // steepsetDescentIter(3,B,r,x);
    // printf("%3d:%3.5f,  %3.5f,  %3.5f\n",iter,x[0],x[1],x[2]);
    // do{
    //     conjugateGradientIter(3,B,r,rTr_iter,d,x);
    //     printf("%3d:%3.5f,  %3.5f,  %3.5f\n",iter, x[0],x[1],x[2]);
    //     norm = inftyNorm(3,r);
    //     iter++;
    // }while(iter<MAXITER && norm/norm0>EPS8);
    // printf("%d\n",iter-1);
    

    for(int k=0;k<len;k++){
        n = nlist[k];
        printf("----- n=%3d -----\n",n);
        A = (double*)malloc(sizeof(double)*n*n);
        f = (double*)malloc(sizeof(double)*n);
        u = (double*)malloc(sizeof(double)*n);
        generateA(n,A);
        generatef(n,f);
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         printf("%3.5f    ",A[i*n+j]);
        //     }
        //     printf("\n");
        // }
        // for(int i=0;i<n;i++){
        //     printf("%3.5f    ",f[i]);
        // }
        // printf("\n");
        memset(u, 0, sizeof(double)*n);

        r = (double*)malloc(sizeof(double)*n);
        d = (double*)malloc(sizeof(double)*n);

        // SD
        for(int i=0;i<n;i++){
            r[i] = f[i];
        }
        norm0 = inftyNorm(n,r);
        iter=1;
        do{
            steepestDescentIter(n,A,r,u);
            norm = inftyNorm(n,r);
            iter++;
            //printf("%3.5f,  %3.5f,  %3.5f\n",u[0],u[1],u[2]);
        }while(iter<MAXITER && norm/norm0>EPS8);
        
        //printf("%3.5f,  %3.5f,  %3.5f\n",u[0],u[1],u[2]);
        printf(">>>SD: %d\n",iter-1);
        //printf("%3.10f\n",norm);


        // CG
        rTr_iter = 0.0;
        memset(u, 0, sizeof(double)*n);  // x0
        for(int i=0;i<n;i++){
            r[i] = f[i];                 // r0
            rTr_iter += r[i]*r[i];       // r0Tr0
            d[i] = r[i];                 // d1
        }
        norm0 = inftyNorm(n,r);
        iter=1;
        steepestDescentIter(n,A,r,u);     // r1,u1
        //printf("%3.5f,  %3.5f,  %3.5f\n",u[0],u[1],u[2]);
        do{
            conjugateGradientIter(n,A,r,&rTr_iter,d,u);
            norm = inftyNorm(n,r);
            iter++;
            //printf("%3.5f,  %3.5f,  %3.5f\n",u[0],u[1],u[2]);
        }while(iter<MAXITER && norm/norm0>EPS8);
        //printf("%3.5f,  %3.5f,  %3.5f\n",u[0],u[1],u[2]);
        printf(">>>CG: %d\n",iter-1);
        
        //printf("%3.10f\n",norm);


        free(A);
        free(f);
        free(u);
        free(r);
        free(d);
    }
}



void generateA(int n, double* A) {
    double Factor = 1.0/(n+1);
    Factor = Factor*Factor*Factor;
	for (int i = 0; i < n; i++) {
		A[i*n + i] = (double)(i + 1)*(n - i)*Factor;
		for (int j = i+1; j < n; j++) {
			A[i*n + j] = (double)(i + 1)*(n - j)*Factor;
			A[j*n + i] = A[i*n + j];
		}
	}
}

double fun(double x){
    return x*(1-x)*exp(x);
}

void generatef(int n, double* f){
    double x;
    for(int i=0;i<n;i++){
        x = (double)(i+1)/(n+1);
        f[i] = fun(x);
    }
}


void steepestDescentIter(int n, double* A, double* r, double* x){
    double alpha;
    double u=0.0,v=0.0;
    double* Ad = (double*)malloc(sizeof(double)*n);
    memset(Ad, 0, sizeof(double)*n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            Ad[i] += A[i*n+j]*r[j];
        }
        u += r[i]*r[i];
        v += r[i]*Ad[i];
    }
    alpha = u/v;
    for(int i=0;i<n;i++){
        x[i] += alpha * r[i];
        r[i] -= alpha * Ad[i];
    }
    free(Ad);
}

// Must use double* rTr_iter
void conjugateGradientIter(int n,double* A, double* r, double* rTr_iter, double* d, double* x){
    double beta,alpha=0.0,rTr=0.0,u;
    double* Ad = (double*)malloc(sizeof(double)*n);
    memset(Ad, 0, sizeof(double)*n);
    for(int i=0;i<n;i++){
        rTr += r[i]*r[i];
    }
    beta = rTr/(*rTr_iter);
    *rTr_iter = rTr;
    for(int i=0;i<n;i++){
        d[i] = r[i] + beta *d[i];
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            Ad[i] += A[i*n+j]*d[j];
        }
        alpha += d[i]*Ad[i];
    }
    alpha = rTr/alpha;
    for(int i=0;i<n;i++){
        x[i] += alpha * d[i];
        r[i] -= alpha * Ad[i];
    }
    free(Ad);
}

double inftyNorm(int n, double *x){
    double res = 0.0;
    for(int i=0;i<n;i++){
        if(fabs(x[i])>res){
            res = fabs(x[i]);
        }
    }
    return res;
}