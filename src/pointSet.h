# pragma once

#include "utility.h"
#include <random>


class PointSet
{
	std::default_random_engine re;          // Use this random engine for all rng related to the point set. This allows to have a fixed seeding per points set.
	std::uniform_real_distribution<> randX; // Uniform random numbers over the whole domain in x direction
	std::uniform_real_distribution<> randY; // Uniform random numbers over the whole domain in x direction

	int n;
	double dHex;                            // Reference spacing of hexagonal packing.
	double rel_dmin = 0.87;                 // Relative minimum distance between points; twice the conflict radius
	double rel_rc = 0.65;                   // Relative maximum coverage radius.
	double sdA = 0.038600518;               // Target standard deviation of cell areas; default is from Schlomer thesis p 64.
	bool allStable = false;                 // To implement termination criteria
	double ONE_X;                           // Make it possible to use another size for toroidal domain.
	double ONE_Y;
	double HALF_X;
	double HALF_Y;
	Point marginBL, marginTR;               // Points within these margins affect points in the main replica
	bool twoTiles = false;                  // true if two tiles are optimized at once

	int maxIterations = 5000;

	struct Replica
	{
		VH vh;
		size_t dtId;
	};

	struct Site
	{
		std::vector<size_t> replicaIds;
		bool isStable = false;
		bool becomeStable = false;
	};

	struct Arrangement
	{
		//size_t dtId;
		std::vector<size_t> replicaIdsToIterate;
	};

	std::vector<Replica> replicas;
	std::vector<Arrangement> arrangements;
	std::vector<Site> sites;
	std::vector<DT> dts;

	Point& getPoint(size_t replicaId);

	/// @return: The point in the center tile
	Point getMainReplica(const Point& p) const;

	void moveSite(size_t siteId, Point targetPoint);
	void moveSite(size_t siteId, Vector shift);

	void coverage(size_t replicaId);
	void conflict(size_t replicaId);
	void capacity(size_t replicaId);

	Point createReplica(Point& p, int i) const;

	/// @return: Randomly ordered vector with the integers from 0 to n-1
	std::vector<size_t> shuffle(const size_t n);

	PointSet(int number_of_points, double aspectRatio);
public:
	// Copying is not allowed since this would invalidate all vertex-, face- and edge-handles when the triangulation is copied
	PointSet(const PointSet&) = delete;
	PointSet operator=(const PointSet&) = delete;

	PointSet(PointSet&&);
	PointSet& operator=(PointSet&&) = default;

	PointSet(int numberOfPoints, int initType, double aspectRatio); // initType:   0: random, 1: darts, 2: jittered grid, 3: regular grid
	PointSet(int nPoints, double* inputPoints, double aspectRatio);
	PointSet(int nPoints, double* inputPoints, double* inputPointsFixedTile, double aspectRatio);

	void setdmin(double d) { rel_dmin = d; };                                   // Set target NND for spring().
	void setRc(double r) { rel_rc = r; };
	void set_sdA(double sd) { sdA = sd; };

	/// @ return: Iteration count until convergence. -1 if no convergence until the specified maximum iteration.
	int ppo();

	void minDistanceCheck();

	void PointSet::getPoints(double* outMatrix)
	{
		if (twoTiles) // points of both point sets are returned in one matrix
		{
			auto arrangement = arrangements[0]; // contains ht points of tile 0 by convention
			auto row = 0;
			for (auto replicaId : arrangement.replicaIdsToIterate)
			{
				Point p = getMainReplica(getPoint(replicaId));
				//printCoordinates(p);
				outMatrix[row] = p.x();
				outMatrix[2 * n + row] = p.y();
				++row;
			}
			arrangement = arrangements[1];
			for (auto replicaId : arrangement.replicaIdsToIterate)
			{
				Point p = getMainReplica(getPoint(replicaId));
				//printCoordinates(p);
				outMatrix[row] = p.x();
				outMatrix[2 * n + row] = p.y();
				++row;
			}
		}
		else
		{
			auto arrangement = arrangements[0]; // contains ht points of tile 0 by convention
			auto row = 0;
			for (auto replicaId : arrangement.replicaIdsToIterate)
			{
				Point p = getMainReplica(getPoint(replicaId));

				outMatrix[row] = p.x();
				outMatrix[n + row] = p.y();
				++row;
			}
		}
	}

	void PointSet::getPoints(Arrangement& arrangement, std::vector<double>& xOut, std::vector<double>& yOut)
	{
		xOut.resize(arrangement.replicaIdsToIterate.size());
		yOut.resize(arrangement.replicaIdsToIterate.size());
		auto i = 0;
		for (auto replicaId : arrangement.replicaIdsToIterate)
		{
			Point p = getMainReplica(getPoint(replicaId));
			xOut[i] = p.x();
			yOut[i] = p.y();
			i++;
		}
	}
};
