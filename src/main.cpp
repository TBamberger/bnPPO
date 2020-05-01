#include "ppo.h"

#include <random>
#include <vector>

int main()
{
	const auto dMin = 0.85;
	const auto rC = 3;// 0.67;
	const auto areaDeltaMax = -1;

	const auto nPoints = 300;
	const auto coordinatesPerPoint = 2;
	std::vector<double> outMatrix(nPoints * coordinatesPerPoint * 2);

	const auto initType = 0;
	const auto aspectRatio = 1.0;
	
	//optimizePattern(dMin, rC, areaDeltaMax, nPoints, initType, outMatrix, aspectRatio);
	
	std::default_random_engine re(12345);
	const auto randX = std::uniform_real_distribution<double>(0, 1);
	std::vector<double> points(nPoints * coordinatesPerPoint);
	for (auto i = 0; i < nPoints * coordinatesPerPoint; ++i)
	{
		points[i] = randX(re);
	}
	//optimizePattern(dMin, rC, areaDeltaMax, nPoints, points, outMatrix, aspectRatio);

	std::vector<double>  points2(nPoints * coordinatesPerPoint);
	for (auto i = 0; i < nPoints * coordinatesPerPoint; ++i)
	{
		points2[i] = randX(re);
	}
	optimizePattern(dMin, rC, areaDeltaMax, nPoints, points.data(), points2.data(), outMatrix.data(), aspectRatio);
}
