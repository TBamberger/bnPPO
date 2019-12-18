#include "ppo.h"

int main()
{
	const auto dMin = 0.85;
	const auto rC = 0.67;
	const auto areaDeltaMax = 0.0;

	const auto nPoints = 300;
	const auto coordinatesPerPoint = 2; // 2D
	auto* outMatrix = new double[nPoints * coordinatesPerPoint];

	const auto initType = 0;
	const auto aspectRatio = 5.0;
	
	optimizePattern(dMin, rC, areaDeltaMax, nPoints, initType, outMatrix, aspectRatio);
	
	delete[] outMatrix;
}
