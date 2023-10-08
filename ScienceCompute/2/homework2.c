#include <stdio.h>
#include <math.h>

#define EPS8 1e-8
#define MAXITER 10000


// Answer to Question 1
int Question1();
double f1(double x);     // Function in Question 1
double Df1(double x);    // The derivative to f1

int main(){
    Question1();
    return 0;
}


double f1(double x){
    return 10.0*x*(1-x)*(x-1.0/4.0)-1.0/4.0;
}
double Df1(double x){
    return (((10.0-12.0*x)*x)-1.0)*5.0/2.0;
}


int Question1(){
    int iter =0;
    double x0=1.0;
    double xk,fk,Dfk;

    printf("========== Question 1 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        fk  = f1(xk);
        Dfk = Df1(xk);
        xk  = xk-fk/Dfk;
        printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);


    printf("\nQustion (a)\n");
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        fk  = f1(xk);
        Dfk = Df1(xk);
        xk  = xk-fk/Dfk/2.0;
        printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);






    return 0;
}