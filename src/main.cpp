#include "ppo.h"

#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

int run(const unsigned seed)
{
	const auto dMin = 0.85;
	const auto rC = 3; // 0.67;
	const auto areaDeltaMax = -1;

	const auto nPoints = 300;
	const auto coordinatesPerPoint = 2;
	std::vector<double> outMatrix(nPoints * coordinatesPerPoint * 2);

	const auto initType = 0;
	const auto aspectRatio = 1.0;

	//optimizePattern(dMin, rC, areaDeltaMax, nPoints, initType, outMatrix, aspectRatio);

	std::default_random_engine re(seed);
	const auto randX = std::uniform_real_distribution<double>(0, 1);
	std::vector<double> points(nPoints * coordinatesPerPoint);
	for (auto i = 0; i < nPoints * coordinatesPerPoint; ++i)
	{
		points[i] = randX(re);
	}
	//optimizePattern(dMin, rC, areaDeltaMax, nPoints, points, outMatrix, aspectRatio);

	std::vector<double> points2(nPoints * coordinatesPerPoint);
	for (auto i = 0; i < nPoints * coordinatesPerPoint; ++i)
	{
		points2[i] = randX(re);
	}
	
	return optimizePattern(dMin, rC, areaDeltaMax, nPoints, points.data(), points2.data(), outMatrix.data(), aspectRatio);
}

template<typename T>
void printMinMeanMax(std::vector<T>& v, const std::string& valueName)
{
	if (v.empty())
	{
		std::cerr << "Warning: Cant' determine min max and mean of empty vector" << std::endl;
		return;
	}
	
	auto mean = std::accumulate(v.begin(), v.end(), 0) / v.size();
	auto [min, max] = std::minmax_element(v.begin(), v.end());

	std::cout << valueName << " (min): " << *min << '\n';
	std::cout << valueName << " (mean): " << mean << '\n';
	std::cout << valueName << " (max): " << *max << std::endl;
}

void testStatistics(const int n)
{
	std::vector<long long> timesMilliseconds(n);
	
	std::vector<int> iterations(n);
	
	for (auto i = 0; i < n; ++i)
	{
		auto startTime = std::chrono::steady_clock::now();
		iterations[i] = run(i);
		auto endTime = std::chrono::steady_clock::now();
		timesMilliseconds[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	}

	printMinMeanMax(iterations, "Iterations");
	printMinMeanMax(timesMilliseconds, "Time [ms]");
}

int main()
{
	testStatistics(100);
}
