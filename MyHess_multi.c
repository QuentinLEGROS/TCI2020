#include "mex.h"
#include <math.h>
#include <limits.h>
//inputs positions 
#define Scaln 0
#define Pn 1
#define Hn 2

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Variables definitions
    //#define coeffout plhs[0]
    int i,j,k;
    int p1,p2;      //P layer indice
    int N,Nx,Nf;      //input data size
    double *X;   //space for interm calculations
    double *P,*Scal,*H,*H_out;
    double temp;
    const size_t *Dim;
    
    // Retrieve data size parameters
    Dim = mxGetDimensions (prhs[Pn]);
    N=Dim[0]; 
    Nx = Dim[1];
    Nf = Dim[2];
//--------- get input parameters --------    
    Scal = mxGetPr(prhs[Scaln]);        //set noisy image pointer
    P = mxGetPr(prhs[Pn]);
    H = mxGetPr(prhs[Hn]);
    //coeffout = mxCreateDoubleMatrix(10,1, mxREAL);
    plhs[0] = mxCreateDoubleMatrix(Nf,Nf, mxREAL);
    H_out = mxGetPr(plhs[0]);
    
//---------- aux variables ---------------

    
    for(i=0;i<Nf*Nf;i++){
        p1 = i%Nf;
        p2 = floor(i/Nf);
        for(j=0;j<N*Nx;j++){
            H[i] -= Scal[j]*P[j+p1*N*Nx]*P[j+p2*N*Nx];
        }
    }
    unsigned int s;
    for(s=0;s<Nf*Nf;s++)
        H_out[s]=H[s];
}

