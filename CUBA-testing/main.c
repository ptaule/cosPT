#include <stdio.h>
#include <math.h>
#include <cuba.h>


static int Integrand(const int *ndim, const cubareal xx[],
        const int *ncomp, cubareal ff[], void *userdata) {

    int a = 10;
    int b = 20;

    ff[0] = sin( a + (b - a)*xx[0] )*cos( a + (b-a)*xx[1] )*exp(1) * pow(b-a,3);

    return 0;
}

#define NDIM 3
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.


int main () {
    int comp, nregions, neval, fail;
    cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

    Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
            EPSREL, EPSABS, VERBOSE | LAST, SEED,
            MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);
    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
            nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                (double)integral[comp], (double)error[comp], (double)prob[comp]);
    return 0;
}

