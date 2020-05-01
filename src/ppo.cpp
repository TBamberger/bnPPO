#include "pointSet.h"

int optimizePattern(double dMin, double rC, double areaDeltaMax, PointSet&& ps, double* outMatrix, double aspectRatio)
{
	ps.setdmin(dMin);
	ps.setRc(rC);
	if (areaDeltaMax >= 0)
		ps.set_sdA(areaDeltaMax);

	const auto nIterations = ps.ppo();
	ps.getPoints(outMatrix);
	return nIterations;
}

int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, int initType, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, initType, aspectRatio);
	return optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}

int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, inMatrix, aspectRatio);
	return optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}

int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* inMatrixFixedTile, double* outMatrix, double aspectRatio)
{
	PointSet ps(nPoints, inMatrix, inMatrixFixedTile, aspectRatio);
	return optimizePattern(dMin, rC, areaDeltaMax, std::move(ps), outMatrix, aspectRatio);
}
