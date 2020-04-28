#include "pointSet.h"

#define FIXED_SEED // If this is defined a fixed seed is used for RNG to achieve reproducible point sets; uses random seed otherwise

// ensure that points don't move out of their tiles!!! Can be accompished by checking if they would moveInTile out of tht tile in the move function
// ensure that I keep the correct points in the center tile!!!

void optimizePattern(double dMin, double rC, double areaDeltaMax, PointSet&& ps, double* outMatrix, double aspectRatio)
{
	ps.setdmin(dMin);
	ps.setRc(rC);
	if (areaDeltaMax >= 0)
		ps.set_sdA(areaDeltaMax);

	ps.ppo();
	ps.getPoints(outMatrix);
}

void optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, int initType, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, initType, aspectRatio);
	optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}

void optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, inMatrix, aspectRatio);
	optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}

void optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* inMatrixFixedTile, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, inMatrix, inMatrixFixedTile, aspectRatio);
	optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}