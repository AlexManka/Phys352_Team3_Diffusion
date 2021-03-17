#include <stdio.h>
#include <math.h>
#include "random.h"

#define step 2.5 
#define MAX 1000
#define walkers 100

int main () {
    int i, j;
    for (i=0; i<walkers; i++) {
        double x[MAX];
        double y[MAX];
        x[0] = 0.0;
        y[0] = 0.0;
        long n = -(i+1);
        for (j=0; j<MAX; j++) {
            double angle = ran1(&n)*2*M_PI;
            x[j+1] = x[j] + step * cos(angle);
            y[j+1] = y[j] + step * sin(angle);
            printf("%f\t%f\t%d\n", x[j], y[j], j);
        }
    }
    return 0;
}