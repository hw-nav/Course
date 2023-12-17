#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288

void Question3();
void Question4();
void Question5();
void Question6();
void Question8();
void Question9();


// g++ -lm main.c
int main(){
    Question3();
    Question4();
    Question5();
    Question6();
    Question8();
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
    double trueValue = 0.8556243918921488031733046202800450612264142850914972603202342815;
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
        printf(">>>%f, %e\n",sum, sum-trueValue);
    }
}



void Question5(){
    printf("================================================\n");
    printf("                   Question 5                   \n");
    printf("================================================\n");
    int n=1, k=5;
    double sum,t;
    double trueValue = 16.0/15.0;
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
        printf(">>>%f, %e\n",sum,sum-trueValue);
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
        }
        sum = 2.0*PI*sum/n;
        printf(">>>%f, %e\n",sum,sum-PI);
    }
}

void Question8(){
    printf("================================================\n");
    printf("                   Question 8                   \n");
    printf("================================================\n");
    int n=1, k=5;
    double trueValue = 1.0-exp(-1.0);
    double sum,t,xc,s1,s2,s3; // xc - x-center
    printf("    |  Simpson  |  Two-Node GL  |  Three-Node GL\n");
    for(int i=0;i<k;i++){
        n = 2*n;
        printf("n=%2d:",n);
        sum = 1.0+exp(-1.0);
        s1=0.0, s2=0.0;
        for(int j=1;j<n;j++){
            s1 += exp(-((double)j-0.5)/n);
            s2 += exp(-(double)j/n);
        }
        s1 += exp(-((double)n-0.5)/n);
        sum += 4.0*s1 + 2.0*s2;
        sum = sum/(6.0*n);
        printf("%f  %e",sum,sum-trueValue);

        sum = 0.0;
        t = sqrt(1.0/3.0)*0.5/n;
        for(int j=1;j<=n;j++){
            xc = (j-0.5)/n;
            sum += exp(-(xc-t)) + exp(-(xc+t));
        }
        sum = sum/(2.0*n);
        printf("   %f  %e",sum,sum-trueValue);

        sum = 0.0;
        t = sqrt(0.6)*0.5/n;
        s1=0.0, s2=0.0, s3=0.0;
        for(int j=1;j<=n;j++){
            xc = (j-0.5)/n;
            s1 += exp(-(xc-t));
            s2 += exp(-xc);
            s3 += exp(-(xc+t));
        }
        sum = 5.0*s1+8.0*s2+5.0*s3;
        sum = sum/(18.0*n);
        printf("   %f  %e\n",sum,sum-trueValue);
    }
}


