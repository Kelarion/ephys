// closest.c
// C = closest(A,B)
// Inputs:  A = full real double vector (reference), sorted ascending
//          B = full real double vector (test),      sorted ascending
// Output:  C = index of closest value in A to value in B
//          i.e., A(C(i)) is the closest value in A to B(i)
//          Ties are resolved in favor of higher index
//          Is not currently coded to handle Inf's and NaN's consistently
//          C is the same size as B
// Programmer:  James Tursa
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t i, j, k, m, n;
    double *A, *B, *C;
//    
// Check input arguments
    if( nrhs != 2 || 
        !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
        mxIsSparse(prhs[0])  || mxIsSparse(prhs[1]) ) {
        mexErrMsgTxt("Need exactly two full double vectors as input");
    }
// Check number of output arguments
    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many outputs");
    }
// Get the number of elements involved
    m = mxGetNumberOfElements(prhs[0]);
    n = mxGetNumberOfElements(prhs[1]);
// Disallow empty reference vector
    if( m == 0 ) {
        mexErrMsgTxt("Reference vector (1st input) cannot be empty");
    }
// Create uninitialized output
    if( mxGetM(prhs[1]) == 1 ) {
        plhs[0] = mxCreateUninitNumericMatrix( 1, n, mxDOUBLE_CLASS, mxREAL );
    } else {
        plhs[0] = mxCreateUninitNumericMatrix( n, 1, mxDOUBLE_CLASS, mxREAL );
    }
// If B is empty, simply return empty result
    if( n == 0 ) {
        return;
    }
// Get the data pointers
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(plhs[0]);
    k = 0;
// Assign 1st index to all values B that are less than 1st A value
    while( k < n && B[k] < *A ) {
        C[k++] = 1.0;
    }
// Step through until B is between two A values, then test for result
    i = 0;
    for( j=k; j<n; j++ ) {
        while( i+1 < m && B[j] >= A[i+1] ) i++;
        if( i+1 == m ) break;
        if( B[j] - A[i] < A[i+1] - B[j] ) {
            C[j] = i + 1;
        } else {
            C[j] = i + 2;
        }
    }
// Assign last index to all values B that are more than last A value
    while( j < n ) {
        C[j++] = m;
    }
}