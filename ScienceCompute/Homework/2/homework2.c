#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>



#define EPS8       1e-8
#define MAXITER    10000
#define MAX_MATRIX 800




// Common Function
void printMatrix(double *A, int row, int col);
int GaussEliminationInverse(int n, double* A,double* invA);
void computeMatrixProduct(double* in1,double* in2, int n1, int n2, int n3, double* out);


// Answer to Problem 1
void Problem1();
double f1(double x);     // Function in Problem 1
double Df1(double x);    // The derivative to f1
double phi1(double x);



// Answer to Problem 2
void Problem2();
double f2(double x);     // Function in Problem 2
double Df2(double x);    // The derivative to f2
double phi2(double x);



// Answer to Problem 3
void Problem3();
double f3(double x);     // Function in Problem 3
double Df3(double x);    // The derivative to f3
double phi3(double x);



// Answer to Problem 4
void Problem4();
double f4(int n, double* z, double x);     // Function in Problem 4
int cmp4(const void *a, const void *b);    // Compare Function for qsort
float randfloat(float a, float b);         // Random Function



// Answer to Probelm 6
int generate6(int n, double* A);                  // Generate Matrix in 1-demension array
void Problem6();



// Answer to Problem 7
double f7(double x);
int generate7(int n, double* A);
void Problem7();


// Answer to Problem 10
int generate10(int  n,double *A);
void Problem10();




int main(){
    Problem1();
    printf("\n\n\n");
    Problem2();
    printf("\n\n\n");
    Problem3();
    printf("\n\n\n");
    Problem4();
    printf("\n\n\n");
    // Problem6();        // Time too long
    printf("\n\n\n");
    // Problem7();        // Undefind reference to exp, need to "gcc homework2.c -0 a.out -lm"
    printf("\n\n\n");
    Problem10();
    return 0;
}


// Common Function
void printMatrix(double *A, int row, int col){
    for(unsigned int i=0;i<row;i++){
        for(unsigned int j=0;j<col;j++){
            if(fabs(A[i*col+j])<EPS8){
                printf("          ");
            }else{
                printf("%5.6f   ",A[i*col+j]);
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







// Functions in Problem 1
double f1(double x){
    return 10.0*x*(1-x)*(x-1.0/4.0)-1.0/4.0;
}
double Df1(double x){
    return (((10.0-12.0*x)*x)-1.0)*2.5;
}
double phi1(double x){
    double gamma=0.5;
    return x-gamma * f1(x) / Df1(x);
}
void Problem1(){
    int iter =1;
    double x0=0.0;
    double xk,fk,Dfk;
    double y0,y1,y2;
    double x[MAXITER];       // Save for compute quadratic or linear convergence
    double e1,e2,e3;

    printf("========== Problem 1 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df1(xk);
        xk  = xk-fk/Dfk;
        fk  = f1(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi1(xk);
        fk  = f1(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }



    printf("\nQustion (c)\n");
    iter = 1;
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        y0 = xk;
        y1 = phi1(y0);
        y2 = phi1(y1);
        xk = (y0*y2-y1*y1)/(y0+y2-2.0*y1);
        fk = f1(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }
        
    
    
}




// Functions in Problem 2
double f2(double x){
    return x*(x*x-1.0)-3.0;
}
double Df2(double x){
    return 3.0*x*x-1.0;
}
double phi2(double x){
    double gamma=0.5;
    return x-gamma * f2(x) / Df2(x);
}
void Problem2(){
    int iter =1;
    double x0=0.0;
    double xk,fk,Dfk;
    double y0,y1,y2;
    double e1,e2,e3,x[MAXITER];

    printf("========== Problem 2 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df2(xk);
        xk  = xk-fk/Dfk;
        fk  = f2(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi2(xk);
        fk  = f2(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        x[iter]=xk;
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }



    printf("\nQustion (c)\n");
    iter = 1;
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        y0 = xk;
        y1 = phi2(y0);
        y2 = phi2(y1);
        xk = (y0*y2-y1*y1)/(y0+y2-2.0*y1);
        x[iter]=xk;
        fk = f2(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }
}




// Functions in Problem 3
double f3(double x){
    return x*x*(x*x-4.0)+4.0;
}
double Df3(double x){
    return 4.0*x*(x*x-2.0);
}
double phi3(double x){
    double gamma=0.5;
    return x-gamma * f3(x) / Df3(x);
}
void Problem3(){
    int iter =1;
    double x0=1.0;
    double xk,fk,Dfk;
    double y0,y1,y2;
    double e1,e2,e3,x[MAXITER];

    printf("========== Problem 3 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f3(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df3(xk);
        xk  = xk-fk/Dfk;
        fk  = f3(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f3(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi3(xk);
        fk  = f3(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }



    printf("\nQustion (c)\n");
    iter = 1;
    xk = x0;
    fk = f3(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        y0 = xk;
        y1 = phi3(y0);
        y2 = phi3(y1);
        xk = (y0*y2-y1*y1)/(y0+y2-2.0*y1);
        fk = f3(xk);
        x[iter]=xk;
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else{
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);
        e1=xk-x[iter-2];
        e2=xk-x[iter-3];
        e3=xk-x[iter-4];
        if(fabs((e1*e3)/(e2*e2)-1.0)<0.5){
            printf("Linear Convergence\n");
        }else if(fabs((e1*e3*e3)/(e2*e2*e2)-1.0)<0.5){
            printf("Quadratic Convergence\n");
        }else{
            printf("Faster than Quadratic Convergence\n");
        }
    }
}




// Functions in Problem 4
double f4(int n, double* z, double x){
    double sum = 0.0;
    for(int i=0;i<n+1;i++)
        sum += 1.0/(x-z[i]);
    return sum;
}
int cmp4(const void *a, const void *b){
    return ((*(double *)a-*(double *)b)>0?1:-1);   // 升序
}
void Problem4(){
    double z[7];       // random points
    double a,b,fa,fb;  // interval for bisectino method and the value
    double x,fx;       // The root x, and f(x)
    int len=4/*length of n*/, n[]={3,4,5,6}, iter;

    printf("========== Problem 4 ==========\n");

    printf("Qustion (b)\n");
    for(int k=0;k<len;k++){
        printf(">>>n=%d\n",n[k]);

        // Generate z[i]
        for(int i=0;i<n[k]+1;i++){
            z[i]=randfloat(0.0,1.0);
            // printf("z%d=%1.4f   ",i,z[i]);
        }
        qsort(z,n[k]+1,sizeof(z[0]),cmp4); // Sort from small to big
        for(int i=0;i<n[k]+1;i++){
            printf("z%d=%1.4f   ",i,z[i]); // Print z[i]
        }
        printf("\n");

        // Bisection Method find the root
        for(int i=0;i<n[k];i++){
            a=z[i],b=z[i+1], fa=f4(n[k],z,a), fb=f4(n[k],z,b);
            x=(a+b)/2;
            fx=f4(n[k],z,x);
            iter = 0;
            while(fabs(fx)>EPS8 && iter<MAXITER){
                if(fa*fx<0){
                    b=x;
                    fb=fx;
                }
                else{
                    a=x;
                    fa=fx;
                }
                x=(a+b)/2;
                fx=f4(n[k],z,x);
                iter++;
            }

            // Test for validation
            //if(iter==MAXITER)
            //    printf("k:MAXITERATE, x=%1.9f, f(x)=%2.9f\n",x,fx);
            //else
            //    printf("k:%d, x=%5.9f, f(x)=%2.9f\n",iter,x,fx);

            if(iter==MAXITER)
                printf("DIVERGENCE  ");       // Error
            else
                printf("x%d=%1.4f   ",i+1,x);
        }
        printf("\n");
    }
}
float randfloat(float a, float b){
    srand((unsigned int)time(0)+(unsigned int)rand());
    return a+(b-a)*(float)rand()/RAND_MAX/(b-a);
}



// Functions in Problem 6
int generate6(int n, double* A) {
	for (int i = 0; i < n; i++) {
		A[i*n + i] = (double)(i + 1)*(n - i)/(n+1);
		for (int j = i+1; j < n; j++) {
			A[i*n + j] = (double)(i + 1)*(n - j)/(n+1);
			A[j*n + i] = A[i*n + j];
		}
	}
	return 0;
}
void Problem6(){
    int len=6/*length of n*/, n[]={4,8,100,200,400,800};
    double *A,*invA;
    clock_t start,end;

    printf("========== Problem 6 ==========\n");

    for(int i=0;i<len;i++){
        A=(double*)malloc(sizeof(double)*n[i]*n[i]);
        invA=(double*)malloc(sizeof(double)*n[i]*n[i]);
        generate6(n[i],A);
        start=clock();
        GaussEliminationInverse(n[i],A,invA);
        end=clock();
        printf(">>>n=%3d, Time=%7.3fs,\n",n[i],(double)(end-start)/CLOCKS_PER_SEC);
        if(n[i]==4 || n[i]==8){
            printf("inv(A)=\n");
            printMatrix(invA,n[i],n[i]);
        }
        
        free(A);
        free(invA);
    } 
}





// Functions in Problem 7
// exp error: Undefined reference to exp
// Need to gcc -o a.out homework2.c -lm
double f7(double x){
    // return x*(x+3)*exp(x);
    return x*(x+3);
}
int generate7(int n,double *A){
    memset(A, 0, sizeof(double)*n*n);
    double factor = (double)(n+1)*(n+1);
    double factor2 = factor*2;
    factor = -factor;
    A[0] = factor2;
    for(int i=0;i<n-1;i++){
        A[i*n+i+1]=factor;
        A[(i+1)*n+i]=factor;
        A[(i+1)*n+(i+1)]=factor2;
    }
}
void Problem7(){
    int len=5, n[]={9,99,199,399,799};
    double *A, *invA, *b, *x, *u;
    FILE* file;
    char filename[20];

    printf("========== Problem 7 ==========\n");
 
    // Gauss Elimination Solve Linear System
    for(int i=0;i<len;i++){
        A   =(double*)malloc(sizeof(double)*n[i]*n[i]);
        invA=(double*)malloc(sizeof(double)*n[i]*n[i]);
        b   =(double*)malloc(sizeof(double)*n[i]);
        x   =(double*)malloc(sizeof(double)*n[i]);
        u   =(double*)malloc(sizeof(double)*n[i]);

        generate7(n[i],A);
        double y=0.1;
        for(int k=0;k<n[i];k++){
            x[k]=(double)(k+1)/(n[i]+1);  // k from 0 to n-1
            b[k]=f7(x[k]);
        }
            
        GaussEliminationInverse(n[i],A,invA);
        computeMatrixProduct(invA,b,n[i],n[i],1,u);
        
        if(n[i]==9){
            printf(">>>n=9, Solution u^T=\n");
            printMatrix(u,n[i],1);
        }

        
        sprintf(filename, "Problem_%d.txt", n[i]);
        file = fopen(filename,"w");
        for(int k=0;k<n[i];k++)
            fprintf(file,"%7.6f,%7.6f\n",x[k],u[k]);
        fclose(file);

        

        free(A);
        free(invA);
        free(b);
        free(x);
        free(u);
    }
}




// Functions in Problem 10
void Problem10(){
    printf("========== Problem10 ==========\n");
    double *A, *b, *true_x, *solve_x,*invA;
    int n = 4, iter=1,flag=0;
    double maximum_norm,t;

    while(iter<MAXITER){
        A   =(double*)malloc(sizeof(double)*n*n);
        invA =(double*)malloc(sizeof(double)*n*n);
        true_x=(double*)malloc(sizeof(double)*n);
        solve_x=(double*)malloc(sizeof(double)*n);
        b=(double*)malloc(sizeof(double)*n);

        generate10(n,A);
        
        maximum_norm=0.0;
        for(int i=0;i<n;i++){
            true_x[i]=1.0;
        }
        computeMatrixProduct(A,true_x,n,n,1,b);
        
        flag=GaussEliminationInverse(n,A,invA);
        if(flag==1){
            break;
        }
        computeMatrixProduct(invA,b,n,n,1,solve_x);
        for(int i=0;i<n;i++){
            t=fabs(solve_x[i]-1.0);
            if(maximum_norm<t)
                maximum_norm=t;
        }

        printf(">>>n=%3d,  maximum_norm=%3.9f\n",n,maximum_norm);
            
        free(A);
        free(invA);
        free(true_x);
        free(solve_x);
        free(b);
        n++;
        iter++;
    }
}

int generate10(int n,double *A){
    double factor;
    for(int i=0;i<n;i++){
        factor = 1.0/(i+1);
        for(int j=0;j<=i;j++){
            A[j*n+(i-j)]=factor;
        }
    }
    for(int i=n;i<2*n-1;i++){
        factor = 1.0/(i+1);
        for(int j=i-n+1;j<n;j++){
            A[j*n+(i-j)]=factor;
        }
    }
    return 0;
}