#include "ppo.h"

#include <random>

int main()
{
	const auto dMin = 0.85;
	const auto rC = 0.67;
	const auto areaDeltaMax = 0.0;

	const auto nPoints = 50;
	const auto coordinatesPerPoint = 2; // 2D
	auto* outMatrix = new double[nPoints * coordinatesPerPoint];

	const auto initType = 0;
	const auto aspectRatio = 1.0;
	
	//optimizePattern(dMin, rC, areaDeltaMax, nPoints, initType, outMatrix, aspectRatio);
	
	std::default_random_engine re(12345);
	const auto randX = std::uniform_real_distribution<double>(0, 1);
	auto* points = new double[nPoints * coordinatesPerPoint];
	for (auto i = 0; i < nPoints * coordinatesPerPoint; ++i)
	{
		points[i] = randX(re);
	}
	optimizePattern(dMin, rC, areaDeltaMax, nPoints, points, outMatrix, aspectRatio);
	
	delete[] outMatrix;
}
