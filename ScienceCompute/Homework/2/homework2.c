#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS8 1e-8
#define MAXITER 10000






// Answer to Problem 1
int Problem1();
double f1(double x);     // Function in Problem 1
double Df1(double x);    // The derivative to f1
double phi1(double x);



// Answer to Problem 2
int Problem2();
double f2(double x);     // Function in Problem 2
double Df2(double x);    // The derivative to f2
double phi2(double x);



// Answer to Problem 3
int Problem3();
double f3(double x);     // Function in Problem 3
double Df3(double x);    // The derivative to f3
double phi3(double x);



// Answer to Problem 4
int Problem4();
double f4(int n, double* z, double x);     // Function in Problem 4


// Random Function
float randfloat(void);
float randfloat(float a, float b);




int main(){
    Problem1();
    printf("\n\n\n");
    Problem2();
    printf("\n\n\n");
    Problem3();
    printf("\n\n\n");
    Problem4();
    printf("%f",randfloat(1.0,2.0));
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
int Problem1(){
    int iter =1;
    double x0=0.0;
    double xk,fk,Dfk;
    double y0,y1,y2;

    printf("========== Problem 1 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df1(xk);
        xk  = xk-fk/Dfk;
        fk  = f1(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f1(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi1(xk);
        fk  = f1(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);



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
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);

    return 0;
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
int Problem2(){
    int iter =1;
    double x0=0.0;
    double xk,fk,Dfk;
    double y0,y1,y2;

    printf("========== Problem 2 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df2(xk);
        xk  = xk-fk/Dfk;
        fk  = f2(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi2(xk);
        fk  = f2(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);



    printf("\nQustion (c)\n");
    iter = 1;
    xk = x0;
    fk = f2(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        y0 = xk;
        y1 = phi2(y0);
        y2 = phi2(y1);
        xk = (y0*y2-y1*y1)/(y0+y2-2.0*y1);
        fk = f2(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);

    return 0;
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
int Problem3(){
    int iter =1;
    double x0=1.0;
    double xk,fk,Dfk;
    double y0,y1,y2;

    printf("========== Problem 3 ==========\n");

    printf("Qustion (a)\n");
    xk = x0;
    fk = f3(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        Dfk = Df3(xk);
        xk  = xk-fk/Dfk;
        fk  = f3(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);


    printf("\nQustion (b)\n");
    iter = 1;
    xk = x0;
    fk = f3(xk);
    while(fabs(fk)>EPS8 && iter<MAXITER){
        xk  = phi3(xk);
        fk  = f3(xk);
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);



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
        // printf("%5d  %5.9f  %5.9f\n",iter,xk,fk);
        iter++;
    }
    if(iter==MAXITER)
        printf("k:MAXITERATE, xk=%5.9f, f(xk)=%5.9f\n",xk,fk);
    else
        printf("k:%d, xk=%5.9f, f(xk)=%5.9f\n",iter,xk,fk);

    return 0;
}




// Functions in Problem 4
double f4(int n, double* z, double x){
    double sum = 0.0;
    for(int i=0;i<n;i++)
        sum += 1.0/(x-z[i]);
    return sum;
}
int Problem4(){
    return 0;
}







// Random Function
float randfloat(void){
    float v=(rand(%__INT_MAX__))/(__INT_MAX__-1.0);
    return v;
}
float randfloat(float a, float b){
    float v=a+(b-a)*randflloat();
    return v;
}