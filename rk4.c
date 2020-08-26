#include "rk4.h"
#include<stdio.h>
static void
rkstep(vecfield_t f, real_t t, real_t h, size_t n, real_t y[n])
{
	real_t dy[n],k1[n],k2[n],k3[n],k4[n], ytemp[n],ttemp;  //ttemp[n] to calculate t values and ytemp[n] to calculate y values of the function
   //  include other arrays as needed
 
   if(n == 0) return;
   
   
  
  f(t, n, y,dy);
     
  for(size_t j = 0; j< n; j++){
      k1[j] = dy[j];   
      ytemp[j] = y[j] + (h*0.5)*k1[j];  //computing 'y' for k2
         
  }

 
  ttemp = t + (h*0.5);            //computing 't' for k2
  f(ttemp, n, ytemp, dy); 
  for(size_t j = 0; j< n; j++){
      
      k2[j] = dy[j];      
      ytemp[j] = y[j] + (h*0.5)*k2[j];  //computing 'y' for k3
      
     
                       
 }

  ttemp = t + (h*0.5);            //computing 't' for k3 
  f(ttemp, n, ytemp, dy);
  for(size_t j = 0; j< n; j++){
      
      k3[j] = dy[j]; 
      ytemp[j] = y[j] + (h)*k3[j];     //computing 'y' for k4
                                 
 }

  ttemp = t + (h);  //computing 't' for k4
  f(ttemp, n, ytemp, dy);
              
  for(size_t j = 0; j< n; j++){
      
      k4[j] = dy[j];   
      y[j] = y[j]+ (1.0/6.0)*(h)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]); //computing next value of y(y[n+1])
      
      }
     // printf("y value in rk 4=%g\n",y);

 
    // your previous code here
    return;
}

// RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal
// time steps.
void
rk4(vecfield_t f, size_t n, real_t y[n], real_t tstart, real_t tend, 
  size_t nsteps)
{
  	real_t h=2.0*PI/256.0;
  	real_t i=0.0;
  	real_t t=0.0;
	 while(i < nsteps)                   
    {
    rkstep(f, t, h, n, y);
    i=i+1.0;
    t=t+h;                                    
    }
  
    return;
}
