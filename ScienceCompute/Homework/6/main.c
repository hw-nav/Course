#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288

void Question3();
void Question4();
void Question5();
void Question6();


// g++ -lm main.c
int main(){
    Question3();
    Question4();
    Question5();
    Question6();
    return 0;
}


void Question3(){
    printf("================================================\n");
    printf("                   Question 3                   \n");
    printf("================================================\n");
    printf(" k & Different Parameters & Forward Errors & Centered Errors\n");
    double h, x=1.0;
    double fneg, f0 = exp(-x*x), fpos;
    double trueDiff = - 2.0 * exp(-1.0);
    double diffForward, diffCentered;
    for(int k=1;k<16;k++){
        h = pow(10.0,-k);
        fneg = exp(-(x-h)*(x-h));
        fpos = exp(-(x+h)*(x+h));
        diffForward = (fpos-f0)/h;
        diffCentered = (fpos-fneg)/(2.0*h);
        printf("%2d   %3.16f   %15.12f   %15.12f\n",k,h,diffForward-trueDiff,diffCentered-trueDiff);
    }
}


void Question4(){
    printf("================================================\n");
    printf("                   Question 4                   \n");
    printf("================================================\n");
    int n=1, k=5;
    double sum;
    for(int i=0;i<k;i++){
        n = 2*n;
        sum = 0.0;
        sum = 0.5*(1.0+exp(-0.5));
        for(int j=1;j<n;j++){
            sum += exp(-(double)j*j/n/n/2);
            //printf("%f\n",sum);
        }
        sum = sum/n;
        printf(">>>%f\n",sum);
    }
}



void Question5(){
    printf("================================================\n");
    printf("                   Question 5                   \n");
    printf("================================================\n");
    int n=1, k=5;
    double sum,t;
    for(int i=0;i<k;i++){
        n = 2*n;
        sum = 0.0;
        //sum = 0.5*(0.0+0.0);
        for(int j=1;j<n;j++){
            t = 2.0*j/n-1.0;
            t = t*t-1.0;
            sum += t*t;
            //printf("%f\n",sum);
        }
        sum = 2.0*sum/n;
        printf(">>>%f\n",sum);
    }
}



void Question6(){
    printf("================================================\n");
    printf("                   Question 6                   \n");
    printf("================================================\n");
    int a=2,b=1;
    double ab = (double)a*b;
    double a2b2_plus=(double)(a*a+b*b);
    double a2b2_minus=(double)(a*a-b*b);
    double I = PI;   // True Value
    int n=1, k=5;
    double sum;
    for(int i=0;i<k;i++){
        n = 2*n;
        sum = 0.5*a/b;
        for(int j=1;j<n;j++){
            sum += ab/(a2b2_plus-a2b2_minus*cos(2.0*PI*j/n));
            //printf("%f\n",sum);
        }
        sum = 2.0*PI*sum/n;
        printf(">>>%f\n",sum);
    }
}


