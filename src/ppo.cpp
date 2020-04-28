#include "pointSet.h"

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
