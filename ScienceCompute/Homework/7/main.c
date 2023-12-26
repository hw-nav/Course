#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288
#define constN 4
#define EPS8 1e-8
#define MAXITER 1e+4

void Question1();
double f1(double x);
void Question2();
double f2(double t,double u);
void Question3();
double f3(double t,double u);
void Question4();



// g++ -lm main.c
int main(){
    Question1();
    //Question2();
    //Question3();
    Question4();
    return 0;
}



double f1(double u){
    double theta = 0.25; // 1/4
    double res = u*(1.0-u)*(u-theta);
    return res;
}


void Question1(){
    printf("=========== Question 1 ===========\n");
    double N = 5.0, dt;
    double u0=0.3, u1, u2, ut;
    double ul[3];
    for(int i=0;i<4;i++){
        // Delta t
        dt = 1.0/N;

        
        u1 = u0;
        for(int j=0;j<N;j++){
            u1 += dt*f1(u1);
        }
        printf("%f",u1);
        if(i>0){
            printf(",%e", u1-ul[0]);
        }
        ul[0] = u1;

        //Backward Euler
        u1 = u0;
        for(int j=0;j<N;j++){
            ut = u1 + dt * f1(u1);
            u2 = u1 + dt * f1(ut);
            while(fabs(u2-ut)>EPS8){
                ut = u2;
                u2 = u1 + dt * f1(ut);
            }
            u1 = u2;
        }
        printf("   %f",u2);
        if(i>0){
            printf(",%e", u2-ul[1]);
        }
        ul[1] = u2;

        // Trapezoidal Method
        u1 = u0;
        for(int j=0;j<N;j++){
            ut = u1 + dt * f1(u1);
            u2 = u1 + 0.5 * dt * ( f1(u1) + f1(ut));
            while(fabs(u2-ut)>EPS8){
                ut = u2;
                u2 = u1 + 0.5 * dt * ( f1(u1) + f1(ut));
            }
            u1 = u2;
        }
        printf("   %f",u2);
        if(i>0){
            printf(",%e", u2-ul[2]);
        }
        ul[2] = u2;

        // double N: 5->10->20->40
        N = 2.0*N;
        printf("\n");     
    }
}

double f2(double t,double u){
    return u-2.0*t/u;
}

void Question2(){
    printf("\n\n\n=========== Question 2 ===========\n");
    double N=5.0, dt, t;
    double u0 = 1.0, u1, u2;
    double k1,k2,k3,k4;
    double trueValue = sqrt(3.0);
    for(int i=0;i<4;i++){
        // Delta t
        dt = 1.0/N;

        // initial u
        u1 = u0;

        // u1: 0->1
        for(int j=0;j<N;j++){
            t  = dt*j; // Not j+1
            k1 = f2(t,u1);
            k2 = f2(t+0.5*dt,u1+0.5*dt*k1);
            k3 = f2(t+0.5*dt,u1+0.5*dt*k2);
            k4 = f2(t+dt,u1+dt*k3);
            u2 = u1 + dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
            u1 = u2;
            //printf("k1=%f, k2=%f, k3=%f, k4=%f, u=%f\n",k1,k2,k3,k4,u2);
        }

        printf("%f, %e",u2,u2-trueValue);

        // double N: 5->10->20->40
        N = 2.0*N;
        printf("\n");
    }
}


void Question3(){
    printf("\n\n\n=========== Question 3 ===========\n");

    int N = 10, iter;
    double dt,t1,t2;
    double *trueValue;
    double u0=1.0,u1,u2,ut;
    double k1,k2,k3,k4;

    for(int i=0;i<3;i++){
        printf("--- N=%d ---\n",N);
        dt = 1.0/N;
        trueValue = (double*)malloc(sizeof(double)*(N+1));
        trueValue[0] = 1.0;
        for(int j=1;j<=N;j++){
            trueValue[j] = exp(-(double)j/N);
        }
        // Forward Euler
        printf("Forward Euler\n");
        u1 = u0;
        for(int j=0;j<N;j++){     // 0,1,2,...,N-1
            t1  = j*dt;
            t2  = (j+1)*dt;
            u2 = u1 + dt * f3(t1,u1);
            u1 = u2;
            //printf("%f,%f,%e\n",t2,u2,u2-trueValue[j+1]);
        }
        // Backward Euler
        printf("Backward Euler\n");
        u1 = u0;
        for(int j=0;j<N;j++){      // 0,1,2,...,N-1
            t1 = j*dt;
            t2 = (double)(j+1)*dt;  // 1,2,3,...,N
            ut = u1 + dt * f3(t1,u1);
            u2 = u1 + dt * f3(t2,ut);
            iter = 0;
            while(fabs(u2-ut)>EPS8 && iter<MAXITER){
                ut = u2;
                u2 = u1 + dt * f3(t2,ut);
                iter++;
            }
            if(iter == MAXITER){
                printf("Error MAX Iterate ---");
            }
            //printf("%f,%f,%e\n",t2,u2,u2-trueValue[j+1]);
            u1 = u2;
        }
        // Trapezoidal Method
        u1 = u0;
        for(int j=0;j<N;j++){     // 0,1,2,...,N-1
            t1 = j*dt;
            t2 = (j+1)*dt;
            ut = u1 + dt * f3(t1,u1);
            u2 = u1 + 0.5 * dt * ( f3(t1,u1) + f3(t2,ut));
            iter = 0;
            while(fabs(u2-ut)>EPS8 && iter<MAXITER){
                ut = u2;
                u2 = u1 + 0.5 * dt * ( f3(t1,u1) + f3(t2,ut));
                iter++;
            }
            if(iter == MAXITER){
                printf("Error MAX Iterate ---");
            }
            // printf("%f,%f,%e\n",t2,u2,u2-trueValue[j+1]);
            u1 = u2;
        }
        // Runge-Kutta
        printf("Runge-Kutta\n");
        u1 = u0;
        for(int j=0;j<N;j++){
            t1  = dt*j; 
            t2  = dt*(j+1);
            k1 = f3(t1,u1);
            k2 = f3(t1+0.5*dt,u1+0.5*dt*k1);
            k3 = f3(t1+0.5*dt,u1+0.5*dt*k2);
            k4 = f3(t1+dt,u1+dt*k3);
            u2 = u1 + dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
            // printf("%f,%f,%e\n",t2,u2,u2-trueValue[j+1]);
            u1 = u2;
            //printf("k1=%f, k2=%f, k3=%f, k4=%f, u=%f\n",k1,k2,k3,k4,u2);
        }

        N = 2*N; // 10->20->40
        printf("\n\n\n");
    }
}

double f3(double t,double u){
    return -10.0 * u + 9.0 * exp(-t);
}



void Question4(){
    printf("\n\n\n=========== Question 4 ===========\n");

    double dt = 0.05;
    int N = 20; // N*dt=1.0
    double x0=PI/6.0,y0=0.0; // x:theta, y:theta prime
    double x1,y1,x2,y2,xt,yt;
    double t1,t2,k1[2],k2[2],k3[2],k4[2];

    // Forward Euler
    printf("Forward Euler\n");
    x1 = x0;
    y1 = y0;
    printf("%f,%f,%f\n",0.0,x1,y1);
    for(int j=0;j<4*N;j++){
        t2 = (j+1)*dt;
        x2 = x1 + y1 * dt;
        y2 = y1 - 16.0 * sin(x1) *dt;
        x1 = x2;
        y1 = y2;
        printf("%f,%f,%f\n",t2,x2,y2);
    }

    // Backward Euler
    printf("Backward Euler\n");
    x1 = x0;
    y1 = y0;
    printf("%f,%f,%f\n",0.0,x1,y1);
    for(int j=0;j<4*N;j++){
        t2 = (j+1)*dt;
        xt = x1 + y1 * dt;
        yt = y1 - 16.0 * sin(x1) * dt;
        x2 = x1 + yt * dt;
        y2 = y1 - 16.0 * sin(xt) * dt;
        while(fabs(xt-x2)>EPS8 || fabs(yt-y2)>EPS8){
            xt = x2;
            yt = y2;
            x2 = x1 + yt * dt;
            y2 = y1 - 16.0 * sin(xt) * dt;
        }
        x1 = x2;
        y1 = y2;
        printf("%f,%f,%f\n",t2,x2,y2);
    }

    // Trapezoidal Method
    printf("Trapezoidal Method\n");
    x1 = x0;
    y1 = y0;
    printf("%f,%f,%f\n",0.0,x1,y1);
    for(int j=0;j<4*N;j++){
        t2 = (j+1)*dt;
        xt = x1 + y1 * dt;
        yt = y1 - 16.0 * sin(x1) * dt;
        x2 = x1 + 0.5 * (y1+yt) * dt;
        y2 = y1 - 8.0 * (sin(x1) + sin(xt)) * dt;
        while(fabs(xt-x2)>EPS8 || fabs(yt-y2)>EPS8){
            xt = x2;
            yt = y2;
            x2 = x1 + 0.5 * (y1+yt) * dt;
            y2 = y1 - 8.0 * (sin(x1) + sin(xt)) * dt;
        }
        x1 = x2;
        y1 = y2;
        printf("%f,%f,%f\n",t2,x2,y2);
    }

    // Runge-Kutta
    printf("Runge-Kutta\n");
    for(int j=0;j<4*N;j++){
        t1  = dt*j; 
        t2  = dt*(j+1);
        k1[0] = y1;
        k1[1] = -16.0*sin(x1);
        k2[0] = y1+0.5*dt*k1[1];
        k2[1] = -16.0*sin(x1+0.5*dt*k1[0]);
        k3[0] = y1+0.5*dt*k2[1];
        k3[1] = -16.0*sin(x1+0.5*dt*k2[0]);
        k4[0] = y1+dt*k3[1];
        k4[1] = -16.0*sin(x1+dt*k3[0]);
        x2 = x1 + dt*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0])/6.0;
        y2 = y1 + dt*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])/6.0;
        x1 = x2;
        y1 = y2;
        printf("%f,%f,%f\n",t2,x2,y2);
    }
}


