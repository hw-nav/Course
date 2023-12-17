#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288
#define constN 4

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
    Question9();
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
        printf("%2d   %e   %e   %e\n",k,h,diffForward-trueDiff,diffCentered-trueDiff);
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



void Question9(){
    printf("================================================\n");
    printf("                   Question 9                   \n");
    printf("================================================\n");

    double I[constN][constN];
    memset(I,0,sizeof(double)*constN*constN); // I[k][log2(N)-2]
    int n=2;

    for(int i=0;i<constN;i++){
        n = 2*n;
        I[i][0] = 0.0;
        for(int j=0;j<n;j++){
            I[i][0] += 4.0/(1.0+(j+0.5)*(j+0.5)/n/n);
        }
        I[i][0] = I[i][0]/n;
        for(int k=1;k<=i;k++){
            I[i][k] = (pow(4.0,k)*I[i][k-1]-I[i-1][k-1])/(pow(4.0,k)-1.0);
        }
    }
    printf(" n=       4          8         16         32\n");
    n = 2;
    for(int k=0;k<constN;k++){
        n = 2*n;
        printf("k=%d",k);
        for(int i=0;i<k;i++){
            printf("           ");
        }
        for(int i=k;i<constN;i++){
            printf("   %f",I[i][k]);
        }
        
        printf("\n");
    }

    printf(" n=         4              8             16             32\n");
    n = 2;
    for(int k=0;k<constN;k++){
        n = 2*n;
        printf("k=%d",k);
        for(int i=0;i<k;i++){
            printf("               ");
        }
        for(int i=k;i<constN;i++){
            printf("   %e",I[i][k]-PI);
        }
        
        printf("\n");
    }


    // double In0[4], In1[3], In2[2], In3[1];
    // int n=2,k=4;
    // printf("      In0      In1      In2     In3\n");
    // for(int i=0;i<k;i++){
    //     n=2*n;
    //     In0[i] = 0.0;
    //     for(int j=0;j<n;j++){
    //         In0[i] += 4.0/(1.0+(j+0.5)*(j+0.5)/n/n);
    //     }
    //     In0[i] =In0[i]/n;
    //     printf("%2d: %f",n, In0[i]);
    //     if(i>0){
    //         In1[i-1] = (4.0*In0[i] - In0[i-1])/3.0;
    //         printf("   %f",In1[i-1]);
    //     }
    //     if(i>1){
    //         In2[i-2] = (16.0*In1[i-1] - In1[i-2])/15.0;
    //         printf("   %f",In2[i-2]);
    //     }
    //     if(i>2){
    //         In3[i-3] = (64.0*In2[i-2] - In2[i-3])/63.0;
    //         printf("   %f",In3[i-3]);
    //     }
    //     printf("\n");
    // }
    // printf("   4   8   16   32\n");
    // printf("en0:%e  %e  %e  %e\n",In0[0]-PI,In0[1]-PI,In0[2]-PI,In0[3]-PI);
    // printf("en1:              %e  %e  %e\n",In1[0]-PI,In1[1]-PI,In1[2]-PI);
    // printf("en2:                            %e  %e\n",In2[0]-PI,In2[1]-PI);
    // printf("en3:                                          %e\n",In3[0]-PI);
}

