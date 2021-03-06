#include "mex.h"
#include "ppo.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

// Input from MATLAB - received via prhs (pass params in this order):
// 0: dMin (which is linked to r_f / conflict radius / min distance between points - relative to r_max,
//        therefore in range [0,1])
// 1: r_c / coverage radius / max distance between a location in the domain and 
//      the nearest sample - relative to r_max, can be >1 but should stay smaller 
//      for decent blue noise properties
// 2: capacity constraint / max delta between areas of voronoi cells of the samples
//    If this value is < 0, the default value of 0.038600518 is used (cf. T. Schl�mer thesis p.64 for justification)
// 3: init type: random, darts, jittered grid, regular grid, specific pattern (passed in the next argument)
// 4: if init type is 'specific': Nx2 input matrix of existing point pattern
//                                Note that a lot of iterations are necessary if the input points are concentrated in a
//                                small region of the input domain because the first iterations are necessary to spread
//                                them according to dMin.
//    else: number of points to generate
// 5: aspect ratio pattern height (y direction) is always 1, pattern width (x direction) is aspectRatio
//
// Output to MATLAB - stored in plhs:
//  Nx2 matrix where each row is a point of the generated pattern
//
// full signature
// noisePattern = blueNoise(dMin, r_c, areaConstraint, initType, nPoints/inputPattern, aspectRatio)

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// Validate and retrieve parameters
	if (nrhs < 6 || nrhs > 7)
		mexErrMsgIdAndTxt("BN:nrhs", "7 inputs required");

	const bool twoTiles = nrhs == 7;

	if (nlhs < 1 || nlhs > 2)
		mexErrMsgIdAndTxt("BN:nlhs", "One or two outputs required");

	for (int i = 0; i < 2; ++i)
		if (!mxIsScalar(prhs[i]) || !mxIsDouble(prhs[i]))
			mexErrMsgIdAndTxt("BN:notScalar", "r_f, r_c and capacity constraints have to be scalar");
	if (!mxIsChar(prhs[3]))
		mexErrMsgIdAndTxt("BN:notString", "initialization type has to be a string / char array");

	const double rF = mxGetScalar(prhs[0]);
	const double rC = mxGetScalar(prhs[1]);
	const double capacityConstraint = mxGetScalar(prhs[2]);
	const double aspectRatio = mxGetScalar(prhs[5]);

	const size_t initTypeSize = mxGetN(prhs[3]) + 1;
	char* initType = new char[initTypeSize];
	if (mxGetString(prhs[3], initType, initTypeSize))
		mexErrMsgIdAndTxt("BN:initType", "Can't retrieve init type");

	int iterations = -1;
	if (std::strcmp(initType, "specific") == 0)
	{
		if (!mxIsDouble(prhs[4]) || mxGetN(prhs[4]) != 2)
			mexErrMsgIdAndTxt("BN:notDouble", "point coordinates have to be a Nx2 matrix of doubles");

		double* inMatrix = mxGetPr(prhs[4]);
		const size_t nColsIn = mxGetN(prhs[4]);
		const size_t nRowsIn = mxGetM(prhs[4]);

		if (twoTiles)
			plhs[0] = mxCreateDoubleMatrix(static_cast<mwSize>(nRowsIn) * 2, static_cast<mwSize>(nColsIn), mxREAL);
		else
			plhs[0] = mxCreateDoubleMatrix(static_cast<mwSize>(nRowsIn), static_cast<mwSize>(nColsIn), mxREAL);
		double* outMatrix = mxGetPr(plhs[0]);
		
		if (twoTiles)
		{
			if (!mxIsDouble(prhs[6]) || mxGetN(prhs[6]) != 2)
				mexErrMsgIdAndTxt("BN:notDouble", "point coordinates have to be a Nx2 matrix of doubles");

			double* inMatrix2 = mxGetPr(prhs[6]);

			if (mxGetM(prhs[6]) != nRowsIn)
				mexErrMsgIdAndTxt("BN:incompatibleTileSize", "fixed and optimized tile need to have the same number of points");

			try
			{
				iterations = optimizePattern(rF, rC, capacityConstraint, nRowsIn, inMatrix, inMatrix2, outMatrix, aspectRatio);
			}
			catch (std::invalid_argument& e)
			{
				mexErrMsgIdAndTxt("BN:invalidInput", e.what());
			}
		}
		else
		{
			try
			{
				iterations = optimizePattern(rF, rC, capacityConstraint, nRowsIn, inMatrix, outMatrix, aspectRatio);
			}
			catch (std::invalid_argument& e)
			{
				mexErrMsgIdAndTxt("BN:invalidInput", e.what());
			}
		}
		// Generate blue noise
	}
	else
	{
		if (!mxIsScalar(prhs[4]))
			mexErrMsgIdAndTxt("BN:notScalar", "if init type is not 'specific' the number of initial points has to be given as a scalar");
		const int nPoints = static_cast<int>(*mxGetPr(prhs[4]));

		// Allocate output
		const size_t nCols = 2; // Only 2D points are supported
		plhs[0] = mxCreateDoubleMatrix(static_cast<mwSize>(nPoints), static_cast<mwSize>(nCols), mxREAL);
		double* outMatrix = mxGetPr(plhs[0]);

		try
		{
			if (std::strcmp(initType, "random") == 0)
			{
				iterations = optimizePattern(rF, rC, capacityConstraint, nPoints, 0, outMatrix, aspectRatio);
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
		catch (std::invalid_argument& e)
		{
			mexErrMsgIdAndTxt("BN:invalidInput", e.what());
		}
	}
	delete[] initType;
	plhs[1] = mxCreateDoubleScalar(static_cast<double>(iterations));
}
