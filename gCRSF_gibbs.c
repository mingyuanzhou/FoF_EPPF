#include "mex.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
/*#include "cokus.c"
#define RAND_MAX_32 4294967295.0*/
        
mwIndex BinarySearch(double probrnd, double *prob_cumsum, mwSize Ksize) {
    mwIndex k, kstart, kend;
    if (probrnd <=prob_cumsum[0])
        return(0);
    else {
        for (kstart=1, kend=Ksize-1; ; ) {
            if (kstart >= kend) {
                /*//k = kend;*/
                return(kend);
            }
            else {
                k = kstart+ (kend-kstart)/2;
                if (prob_cumsum[k-1]>probrnd && prob_cumsum[k]>probrnd)
                    kend = k-1;
                else if (prob_cumsum[k-1]<probrnd && prob_cumsum[k]<probrnd)
                    kstart = k+1;
                else
                    return(k);
            }
        }
    }
    return(k);
}


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    mwSize m;
    mwIndex i, j, k, Ksize, WordNum;
    double r, a, p, mass;
    double *n_k, *z, *SampleDex;
    
    double *prob_cumsum;
    double cum_sum,probrnd;
    void *newptr;

    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxDuplicateArray(prhs[1]);
    
    z = mxGetPr(plhs[0]);
    n_k = mxGetPr(plhs[1]);
    
    SampleDex = mxGetPr(prhs[2]);
    r = mxGetScalar(prhs[3]);
    a = mxGetScalar(prhs[4]);
    p = mxGetScalar(prhs[5]);
    
    
    Ksize = mxGetM(prhs[1])*mxGetN(prhs[1]);
    WordNum = mxGetM(prhs[2])*mxGetN(prhs[2]);
    prob_cumsum = (double *) mxCalloc(Ksize,sizeof(double));
    
    mass = r*pow(p,-a); 
    for (i=0;i<WordNum;i++){
        j = (mwIndex) SampleDex[i] -1;
        k = (mwIndex) z[j] -1;
        if(z[j]>0){
            n_k[k]--;
        }
        for (cum_sum=0, k=0; k<Ksize; k++) {
            cum_sum += n_k[k]-a;
            prob_cumsum[k] = cum_sum;
        }
        
        if ( ((double) rand() / RAND_MAX * (cum_sum + mass) < cum_sum)){
            probrnd = (double)rand()/(double)RAND_MAX*cum_sum;
            k = BinarySearch(probrnd, prob_cumsum, Ksize);
        }
        else{
            for (k=0; k<Ksize; k++){
                if ((mwIndex)n_k[k]==0){
                    break;
                }
            }
            if (k==Ksize){
                Ksize++;
                newptr = mxRealloc(n_k,sizeof(*n_k)*Ksize);
                mxSetPr(plhs[1], newptr);
                mxSetM(plhs[1], Ksize);
                mxSetN(plhs[1], 1);
                n_k = mxGetPr(plhs[1]);
                n_k[Ksize-1]=0;
                prob_cumsum =  mxRealloc(prob_cumsum,sizeof(*prob_cumsum)*Ksize); 
            }
        }
        z[j] = k+1;
        n_k[k]++;
    }
    mxFree(prob_cumsum);
}
