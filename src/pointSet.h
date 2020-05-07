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
	bool twoTiles = false;                  // true if two tiles are optimized at once

	int maxIterations = 500;

	inline const static double OFFSET_LENGTH = 1e-10; /// If input contains identical points or points are would be moved to the same position, they are moved by this distance to keep points separated.

	struct Replica
	{
		VH vh;
		size_t dtId = -1;
	};

	struct Site
	{
		std::vector<size_t> replicaIds;
		bool isStable = false;
		bool becomeStable = false;
	};

	struct Arrangement
	{
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

	PointSet(int nPoints, double aspectRatio);

	/// @return: Random vector with the specified length generated with the random engine of the point set
	Vector randomVector(double vectorLength);
	VH insertUnique(DT& dt, const Point& p);

public:
	// Copying is not allowed since this would invalidate all vertex-, face- and edge-handles when the triangulation is copied
	PointSet(const PointSet&) = delete;
	PointSet operator=(const PointSet&) = delete;

	PointSet(PointSet&&) = default;
	PointSet& operator=(PointSet&&) = default;

	PointSet(int nPoints, int initType, double aspectRatio); // initType:   0: random, 1: darts, 2: jittered grid, 3: regular grid
	PointSet(int nPoints, double* inputPoints, double aspectRatio);
	PointSet(int nPoints, double* inputPoints, double* inputPoints2, double aspectRatio);

	~PointSet() = default;

	void setDMin(double d); // Set target NND for spring().
	void setRc(double r);
	void setSdA(double sd);

	/// @ return: Iteration count until convergence. -1 if no convergence until the specified maximum iteration.
	int ppo();

	void minDistanceCheck();

	void getPoints(double* outMatrix);

	void getPoints(Arrangement& arrangement, std::vector<double>& xOut, std::vector<double>& yOut);

	[[nodiscard]] bool isInMainReplica(const Point& p) const;
};
