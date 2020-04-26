# pragma once

#include "utility.h"
#include <random>


class PointSet
{
	std::default_random_engine re;          // Use this random engine for all rng related to the point set. This allows to have a fixed seeding per points set.
	std::uniform_real_distribution<> randX; // Uniform random numbers over the whole domain in x direction
	std::uniform_real_distribution<> randY; // Uniform random numbers over the whole domain in x direction

	// todo: unchecked:
	int n;
	double dhex;                            // Reference spacing of hexagonal packing.
	double rel_dmin = 0.87;                 // Relative minimum distance between points; twice the conflict radius
	double rel_rc = 0.65;                   // Relative maximum coverage radius.
	double sdA = 0.038600518;               // Target standard deviation of cell areas; default is from Schlomer thesis p 64.
	bool allStable = false;                 // To implement termination criteria
	double ONE_X;                           // Make it possible to use another size for toroidal domain.
	double ONE_Y;
	double HALF_X;
	double HALF_Y;
	Point marginBL, marginTR;                                                   // Points within these margins affect points in the main replica

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

	// todo: consider using fixed size containers where applicable
	std::vector<Replica> replicas;
	std::vector<Arrangement> arrangements;
	std::vector<Site> sites;
	std::vector<DT> dts;

	Point& getPoint(size_t replicaId);

	/// Returns the point in the center tile
	Point getMainReplica(const Point& p) const;
	
	void moveSite(size_t siteId, Point targetPoint);
	void moveSite(size_t siteId, Vector shift);
	
	void coverage(size_t replicaId);
	void conflict(size_t replicaId);
	void capacity(size_t replicaId);

	//Point PointSet::setSite(int index, Point p) {                        // Set location of the indexed point (in t-domain) and insert it in triangulation
	//	p = getMainReplica(p);
	//	sites[index].p = p;                                               // Save a handy copy of point coordinates
	//	sites[index].isStable = false;
	//	for (int i = 0; i < 9; i++)
	//	{
	//		sites[index].vh[i] = dt.insert(createReplica(p, i));       // insert replica in triangulation and keep handle to it
	//		sites[index].vh[i]->info().id = index;                        // Point the DT point back to map entry
	//	}
	//	return p;
	//}

	Point createReplica(Point& p, int i) {                                    // Find one of the 9 replicas of a point
		i = (i + 4) % 9;                                                // We make the middle replica at index 0
		double x = p.x() + (i % 3 - 1) * ONE_X;                         // Add -ONE, 0, or ONE to x
		double y = p.y() + (i / 3 - 1) * ONE_Y;                         // Same for y
		return Point(x, y);
	}

	PointSet(int number_of_points, double aspectRatio);
public:
	// Copying is not allowed since this would invalidate all vertex-, face- and edge-handles when the triangulation is copied
	PointSet(const PointSet&) = delete;
	PointSet operator=(const PointSet&) = delete;

	PointSet(PointSet&&);
	PointSet& operator=(PointSet&&) = default;

	PointSet(int number_of_points, int initType, double aspectRatio); // initType:   0: random, 1: darts, 2: jittered grid, 3: regular grid
	PointSet(int nPoints, double* inputPoints, double aspectRatio);
	PointSet(int nPoints, double* inputPoints, double* inputPointsFixedTile, double aspectRatio);

	void setdmin(double d) { rel_dmin = d; };                                   // Set target NND for spring().
	void setRc(double r) { rel_rc = r; };
	void set_sdA(double sd) { sdA = sd; };
	
	/// @ return: Iteration count until convergence. -1 if no convergence until the specified maximum iteration.
	int ppo();

	void PointSet::getPoints(double* outMatrix)
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
};