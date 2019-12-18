#include "mex.h"
#include "ppo.h"
#include <cstring>
#include <iostream>


// rhs (pass params in this order):
// 0: dMin (which is linked to r_f / conflict radius / min distance between points - relative to r_max,
//        therefore in range [0,1])
// 1: r_c / coverage radius / max distance between a location in the domain and 
//      the nearest sample - relative to r_max, can be >1 but should stay smaller 
//      for decent blue noise properties
// 2: capacity constraint / max delta between areas of voronoi cells of the samples
// 3: init type: random, darts, jittered grid, regular grid, specific pattern (passed in the next argument)
// 4: if init type is 'specific': Nx2 input matrix of existing point pattern
//    else: number of points to generate
// 5: aspect ratio pattern height (y direction) is always 1, pattern width (x direction) is aspectRatio
//
// lhs:
//  Nx2 matrix where each row is a point of the generated pattern
//
// full signature
// noisePattern = myPPO(dMin, r_c, d_min?, initType, nPoints/inputPattern, aspectRatio

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Validate and retrieve parameters
    if (nrhs != 6)
        mexErrMsgIdAndTxt("BN:nrhs", "6 inputs required");
    
    if (nlhs != 1)
        mexErrMsgIdAndTxt("BN:nlhs", "One output required");
    
    for (int i = 0; i < 2; ++i)
        if (!mxIsScalar(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgIdAndTxt("BN:notScalar", "r_f, r_c and capacity constraints have to be scalar");
    if (!mxIsChar(prhs[3]))
        mexErrMsgIdAndTxt("BN:notString", "initialization type has to be a string / char array");
    
    double rF = mxGetScalar(prhs[0]);
    double rC = mxGetScalar(prhs[1]);
    double capacityConstraint = mxGetScalar(prhs[2]);
    double aspectRatio = mxGetScalar(prhs[5]);
    
    size_t initTypeLen = mxGetN(prhs[3]) + 1;
    char *initType = new char[initTypeLen];
    if (mxGetString(prhs[3], initType, initTypeLen))
        mexErrMsgIdAndTxt("BN:initType", "Can't retrieve init type");
    
    if (std::strcmp(initType, "specific") == 0)
    {
        if (!mxIsDouble(prhs[4]) || mxGetN(prhs[4]) != 2)
            mexErrMsgIdAndTxt("BN:notDouble", "point coordinates have to be a Nx2 matrix of doubles");
    
        double *inMatrix = mxGetPr(prhs[4]);
        size_t nColsIn = mxGetN(prhs[4]);
        size_t nRowsIn = mxGetM(prhs[4]);
        for (mwIndex i = 0; i < nRowsIn; ++i)
        {
            mwIndex subs[2] = {i, 0};
            mwIndex xValIndex = mxCalcSingleSubscript(prhs[4], 2, subs);
            printf("index %d, element (%f %f)\n", xValIndex, inMatrix[xValIndex], inMatrix[xValIndex]);
        }
        
    // Generate blue noise
    }
    else 
    {
        if (!mxIsScalar(prhs[4])) // todo: assert that prhs[4] is an integer value
            mexErrMsgIdAndTxt("BN:notScalar", "if init type is not 'specific' the number of initial points has to be given as a scalar");
        int nPoints = (int)*mxGetPr(prhs[4]);

        // Allocate output
        size_t ncols = 2; // Currently only 2D points are supported
        plhs[0] = mxCreateDoubleMatrix((mwSize)nPoints, (mwSize)ncols, mxREAL);
        double* outMatrix = mxGetPr(plhs[0]);
        
        if (std::strcmp(initType, "random") == 0)
        {
            optimizePattern(rF, rC, capacityConstraint, nPoints, 0, outMatrix, aspectRatio);
        }
        else if (std::strcmp(initType, "darts") == 0)
        {
            mexErrMsgIdAndTxt("BN:initType", "Init type not yet supported");
        }
        else if (std::strcmp(initType, "jittered") == 0)
        {
            mexErrMsgIdAndTxt("BN:initType", "Init type not yet supported");
        }
        else if (std::strcmp(initType, "grid") == 0)
        {
            mexErrMsgIdAndTxt("BN:initType", "Init type not yet supported");
        }
        else
        {
            mexErrMsgIdAndTxt("BN:initType", "Unknown init type");
        }
    }    
}