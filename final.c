#include <stdio.h>
#include <math.h>
#include "random.h"

#define step 2.5 // mean free path of particles in water, in units of Aº
#define MAX 1000 //number of steps taken 
#define walkers 500 //number of particles

/* Main function: initializes the correct number of particles and allows each 
particle to take the correct number of random steps. 
The output prints out MAX * walkers rows of data, with all of the steps for one walker 
followed by all of the steps for the next walker. For example, with MAX = 1000 and walkers = 500, 
rows 1-1000 are for walker #1, rows 1001-2000 are for walker #2, etc. */ 

int main () {
    int i, j;

    /* this outer loop initializes particles and loops over the total number of particles */ 
    for (i=0; i<walkers; i++) {

        /* n1 and n2 are negative integers used to generate random values between 0 and 1 
        for calculating initial position */    
        long n1 = -2*i; 
        long n2 = -2*i-1;

        /* arrays used to store the x and y values of the particle being initialized in each iteration */ 
        double x[MAX];
        double y[MAX];

        /* uses two random integers between 0 and 1 to calculate an initial x and y position 
        for the particle according to a Gaussian distribution (Box-Muller) */
        double radius0 = sqrt(-2*log(ran1(&n1)));
        double angle0 = ran1(&n2)*2*M_PI;
        x[0] = radius0 * cos(angle0);
        y[0] = radius0 * sin(angle0);

        /* n is a negative integer used to generate random values between 0 and 1 
        for calculating the direction of each step */ 
        long n = -(i+1);

        /* this inner loop takes MAX random steps for each particle */ 
        for (j=0; j<MAX; j++) {

            /* random angle between 0 and 2pi, direction of step */ 
            double angle = ran1(&n)*2*M_PI;

            /* calculates next x and y position after a step with fixed step size and direction given by angle */ 
            x[j+1] = x[j] + step * cos(angle);
            y[j+1] = y[j] + step * sin(angle);

            /* prints out the x and y position and the step number for the particle */ 
            printf("%f\t%f\t%d\n", x[j], y[j], j);
        }
    }
    return 0;
}