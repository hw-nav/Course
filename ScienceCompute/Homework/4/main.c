#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define EPS8 1e-8
#define MAXITER 5000000

void generateA(int n, double* A);
double fun(double x);
void generatef(int n, double* f);
void steepsetDescentIter(int n, double* A, double* r, double* x);
void conjugateGradientIter(int n,double* A, double* r, double rTr_iter, double* d, double* x);
double inftyNorm(int n, double *x);


// ERROR: g++ -lm main.c
int main(){
    int len=6, nlist[]={9,19,39,79,159,319};
    int n;
    int iter;
    double *A, *f, *u;
    double *r,*d,rTr_iter;
    double norm, norm0;
    for(int k=0;k<len;k++){
        n = nlist[k];
        A = (double*)malloc(sizeof(double)*n*n);
        f = (double*)malloc(sizeof(double)*n);
        u = (double*)malloc(sizeof(double)*n);
        generateA(n,A);
        generatef(n,f);
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
            steepsetDescentIter(n,A,r,u);
            norm = inftyNorm(n,r);
            iter++;
        }while(iter<MAXITER && norm/norm0>EPS8);
        printf("%d\n",iter);


        // CG
        rTr_iter = 0.0;
        for(int i=0;i<n;i++){
            r[i] = f[i];
            rTr_iter += r[i]*r[i];
            d[i] = r[i];
        }
        norm0 = inftyNorm(n,r);
        iter=1;
        steepsetDescentIter(n,A,r,u);
        do{
            conjugateGradientIter(n,A,r,rTr_iter,d,u);
            norm = inftyNorm(n,r);
            iter++;
        }while(iter<MAXITER && norm/norm0>EPS8);
        printf("%d\n",iter);


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


void steepsetDescentIter(int n, double* A, double* r, double* x){
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


void conjugateGradientIter(int n,double* A, double* r, double rTr_iter, double* d, double* x){
    double beta,alpha=0.0,rTr,u;
    double* Ad = (double*)malloc(sizeof(double)*n);
    memset(Ad, 0, sizeof(double)*n);
    for(int i=0;i<n;i++){
        rTr += r[i]*r[i];
    }
    beta = rTr/rTr_iter;
    rTr_iter = rTr;
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