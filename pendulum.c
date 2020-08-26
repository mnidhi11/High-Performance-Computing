#include "rk4.h"
#include "pendulum.h"
#include <stdlib.h>
#include <stdio.h>

// This code defines the Poincare' map for the forced damped nonlinear
// pendulum, x'' + damping*x' + sin(x) = force*cos(t).

//----------------------------------------------------------------------------
//  The first-order system corresponding to the forced damped pendulum.


static void
pendulum(real_t t, size_t n,  const real_t y[n], real_t dy[n])
{
    for (int i=0;i<n;i++){

        if (i%2==0){

            dy[i]=y[i+1]; //equation 1 applied to even index elements

        }

        else{

            dy[i]=force*cos(t) - sin(y[i-1]) - delta*y[i];  //equation 2 applied to odd index elements

        }
    }

   return;
}

//----------------------------------------------------------------------------
// Poincare map for the forced nonlinear pendulum.  
// Y is the 2-vector giving the angular (position, velocity) of the pendulum after every 2*PI time interval.
void
poincare(int n,real_t y[], size_t niter) // this calculates poincare map 
{
   const size_t NSTEPS = 256;
 
   //printf("value of damping=%g and force=%g\n",delta,force);

   for(size_t j = 0; j < niter; j++)  {

      rk4(pendulum, n, y, 0.0, TWOPI, NSTEPS);

  //wrapping the values of position -pi to pi

       for (int i=0;i<n;i++){

           if ((y[i]<-PI) || (y[i]>PI)) {

            y[i]=remainder(y[i],TWOPI);

           }
       
       }

    } 
 

   return;
}

