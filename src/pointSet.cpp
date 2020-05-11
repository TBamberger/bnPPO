#define FIXED_SEED // If this is defined a fixed seed is used for RNG to achieve reproducible point sets; uses random seed otherwise
//#define PLOT_POINTS // If defined points of arrangement one are plotted each iteration.
//#define PRINT_INFO // If defined additional status information are printed.

#include "pointSet.h"

#include <cassert>
#include <stdexcept>

#ifdef PLOT_POINTS
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

Point& PointSet::getPoint(const size_t replicaId)
{
	return replicas[replicaId].vh->point();
}

Point PointSet::getMainReplica(const Point& p) const
{
	double x = p.x(), y = p.y();
	while (x < 0) x += ONE_X;
	while (x >= ONE_X) x -= ONE_X;
	while (y < 0) y += ONE_Y;
	while (y >= ONE_Y) y -= ONE_Y;
	return Point(x, y);
}

void PointSet::addSite(Point p, std::vector<ReplicaPrototype> replicaPrototypes)
{
#ifndef NDEBUG
	const auto rpWithMaxDtId = std::max_element(replicaPrototypes.begin(), replicaPrototypes.end(), [](const auto& r1, const auto& r2) {return r1.getDtId() < r2.getDtId(); });
	const auto maxDtId = rpWithMaxDtId->getDtId();
	assert(maxDtId < dts.size());
	assert(maxDtId < arrangements.size());
#endif

	p = getMainReplica(p);
	Site site;
	const auto siteId = sites.size();

	// find unique location of p (there is certainly a more efficient method to do this, but this has to do for now due to time constraints)
	if (!replicaPrototypes.empty())
	{
		DT& dt = dts[replicaPrototypes[0].getDtId()];
		const auto initialVertexCount = dt.number_of_vertices();
		dt.insert(p); // In case of an existing point nothing is inserted, therefore no need to remove p again; If p doesn't exist it will be inserted later on anyway
		while (dt.number_of_vertices() == initialVertexCount)
		{
			// If vertex count didn't change, there was an existing point at position p and the vertex handle of the existing
			// point was returned. To keep points separated, the new point is inserted with a tiny offset next to the
			// existing point. This doesn't have a significant impact on the result of the algorithm since points with
			// identical positions will be moved apart from each other to satisfy r_f anyway.
#ifdef PRINT_INFO
			std::cout << "duplicate input points moved away from each other" << std::endl;
#endif
			const Vector offset = randomVector(OFFSET_LENGTH);
			p = getMainReplica(p + offset); // todo: don't get main replica since this relies on the assumption that replicaPrototypes[0] is the main one
			dt.insert(p);
		}
	}
	
	for (auto& rp : replicaPrototypes)
	{
		Replica r;
		const auto replicaId = replicas.size();
		r.vh = dts[rp.getDtId()].insert(createReplica(p, rp.getReplicaNumber()));
		r.vh->info().siteId = siteId;
		r.dtId = rp.getDtId();
		replicas.push_back(r);
		site.replicaIds.push_back(replicaId);
		if (rp.isMainReplica())
		{
			arrangements[rp.getDtId()].replicaIdsToIterate.push_back(replicaId);
		}
	}
	sites.push_back(site);
}

void PointSet::moveSite(size_t siteId, Point targetPoint)
{
	assert(0 <= siteId && siteId < sites.size());
	const auto currentPoint = getMainReplica(getPoint(sites[siteId].replicaIds[0])); // todo: getMainReplica könnte ich mir hier sparen wenn ich sicherstelle, dass siteId 0 immer in der mittleren Kachel ist
	targetPoint = getMainReplica(targetPoint);
	const Vector shift = targetPoint - currentPoint;
	assert(isInMainReplica(currentPoint + shift));

	auto replicaCount = 0;
	for (auto replicaId : sites[siteId].replicaIds)
	{
		auto& replica = replicas[replicaId];
		auto& dt = dts[replica.dtId];
		auto newPosition = replica.vh->point() + shift;

		// new check for existing point
		DT::Locate_type locateType;
		int locateIndex;
		auto fh = dt.locate(newPosition, locateType, locateIndex);
		if (locateType == DT::Locate_type::VERTEX)
		{ // new position coincides with existing vertex
			// There is already an existing point at the target position, if the point would be inserted there the
			// points would be fused together since a delaunay triangulation obviously can't deal with multiple points
			// with identical coordinates.
			// Therefore a minute offset is added to the target position to keep the points separate. This doesn't
			// influence the algorithm performance since points with identical position would be pushed apart anyway to
			// satisfy the conflict radius.
#ifdef PRINT_INFO
			std::cout << "move to existing position avoided" << std::endl;
#endif //PRINT_INFO
			moveSite(siteId, newPosition + randomVector(OFFSET_LENGTH));
			return;
		}
		replica.vh = dt.move(replica.vh, newPosition);

		// Mark neighbors unstable (has to be done for each replica since replicas might be in different triangulations
		// and therefore have different neighbors).
		VC vc = dt.incident_vertices(replica.vh), done(vc);
		do
		{
			if (!dt.is_infinite(vc))
			{
				assert(0 <= vc->info().siteId && vc->info().siteId < sites.size());
				sites[vc->info().siteId].becomeStable = false;
			}
		} while (++vc != done);
		++replicaCount;
	}
	sites[siteId].becomeStable = false;
	allStable = false;
}

void PointSet::moveSite(size_t siteId, const Vector shift)
{
	moveSite(siteId, getPoint(sites[siteId].replicaIds[0]) + shift);
}

void PointSet::coverage(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	const auto rc = rel_rc * dHex;

	const auto vhDegree = replica.vh->degree();
	std::vector<int> id(vhDegree);
	std::vector<double> scale(vhDegree);
	std::vector<Vector> edge(vhDegree);
	auto m = 0;

	FC fc = dt.incident_faces(replica.vh);
	const FC done(fc);
	VC vc = dt.incident_vertices(replica.vh, fc);
	do
	{
		++vc;
		edge[m] = p - vc->point();
		id[m] = vc->info().siteId;
		if (triangleType(fc) <= 0) // acute or right angled
		{
			Point c = dt.circumcenter(fc);
			const double l = length(c - p);
			scale[m] = rc / l;
		}
		else
		{
			scale[m] = 2; // > 1
		}
		m++;
	} while (++fc != done);

	for (auto i = 0u; i < vhDegree; i++)
	{
		const double scl = std::min(scale[i], scale[(i + 1u) % vhDegree]);
		if (scl < 1)
		{
			const Vector shift = (1 - scl) * edge[i];
			moveSite(id[i], shift);
		}
	}
}

void PointSet::conflict(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	const double dMin = rel_dmin * dHex;

	const auto vhDegree = replica.vh->degree();
	std::vector<bool> conflict(vhDegree);
	std::vector<Vector> shift(vhDegree);
	std::vector<int> id(vhDegree);

	int m = 0;
	VC vc = dt.incident_vertices(replica.vh);
	const VC done(vc);
	do
	{
		Vector edge = vc->point() - p;
		const double l = length(edge);
		if (l < dMin)
		{
			conflict[m] = true;
			shift[m] = (1.001 * dMin / l - 1) * edge;
			id[m] = vc->info().siteId;
		}
		else
		{
			conflict[m] = false;
		}
		m++;
	} while (++vc != done);

	for (auto i = 0u; i < vhDegree; i++)
	{
		if (conflict[i]) moveSite(id[i], shift[i]);
	}
}

void PointSet::capacity(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	const auto vhDegree = replica.vh->degree();
	std::vector<double> d(vhDegree);    // Distance to neighbor (= 2 x distance to Voronoi edge)
	std::vector<double> el(vhDegree);   // length of the voronoi edges
	std::vector<Vector> dir(vhDegree);  // Direction vectors to neighbors
	std::vector<int> id(vhDegree);      // Id's of neighbors. We can't use the circulator for updating
	double area = 0; // Area of Voronoi cell
	FC fc2 = dt.incident_faces(replica.vh), fc1(fc2++);               // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
	VC vc = dt.incident_vertices(replica.vh, fc2);                    // The vertex sharing fc1 and fc2 with v[i].vh
	const VC done(vc);
	int m = 0;
	do
	{
		Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);   // Circumcenters of faces are endpoints of Voronoi cell edge
		area += c1.x() * c2.y() - c1.y() * c2.x();                    // Accumulate areas (rhs is cross product)
		el[m] = length(c2 - c1);                                    // Length of Voronoi edge
		dir[m] = (vc->point() - p);
		d[m] = length(dir[m]);
		dir[m] = dir[m] / d[m];                                       // Normalize direction vector
		id[m] = vc->info().siteId;
		++fc1;
		++fc2;
		++m;
	} while (++vc != done);
	area /= 2;
	const double averageArea = ONE_X * ONE_Y / n;
	const double dA = area - averageArea; // Required expansion or contraction
	if (fabs(dA) > sdA)
	{
		double sum_w = 0;
		for (auto j = 0u; j < vhDegree; j++)
		{
			sum_w += el[j] * el[j];
		}
		const double pressure = -2 * dA / sum_w; // pressure per unit length of edges
		for (auto j = 0u; j < vhDegree; j++)
		{
			const Vector force = pressure * el[j] * dir[j];
			moveSite(id[j], force);
		}
	}
}

Point PointSet::createReplica(Point& p, int i) const
{                                                 // Find one of the 9 replicas of a point
	i = (i + 4) % 9;                              // We make the middle replica at index 0
	const double x = p.x() + (i % 3 - 1) * ONE_X; // Add -ONE, 0, or ONE to x
	const double y = p.y() + (i / 3 - 1) * ONE_Y; // Same for y
	return Point(x, y);
}

std::vector<size_t> PointSet::shuffle(const size_t n)
{
	const auto randMax = std::uniform_int_distribution<size_t>();
	std::vector<size_t> v;
	v.resize(n);
	for (size_t i = 0; i < n; ++i)
	{
		v[i] = i;
	}
	for (size_t i = 0; i < n - 1; ++i)
	{
		std::swap(v[i], v[i + randMax(re) % (n - 1 - i)]);
	}
	return v;
}

PointSet::PointSet(const int nPoints, const double aspectRatio) :
#ifdef FIXED_SEED
	re(12345)
#else
	re(std::random_device{}())
#endif // FIXED_SEED
	, n(nPoints)
{
	if (nPoints < 2)
	{
		throw std::invalid_argument("Can't work on less than two points.");
	}
	ONE_X = aspectRatio;
	ONE_Y = 1;
	HALF_X = 0.5 * ONE_X;
	HALF_Y = 0.5 * ONE_Y;
	randX = std::uniform_real_distribution<>(0, ONE_X);
	randY = std::uniform_real_distribution<>(0, ONE_Y);
	dHex = sqrt(ONE_X * ONE_Y * 2 / (sqrt(3) * nPoints)); // Maximum packing distance
}

Vector PointSet::randomVector(const double vectorLength)
{
	// Note that due to the non zero lower bound the random vectors are not perfectly uniformly distributed.
	// This is not required in this application. If it is desired anyhow one could generate in [0,1) and reject the
	// zero vector.
	const static auto randDist = std::uniform_real_distribution<>(0.1, 1);
	Vector v(randDist(re), randDist(re));
	return v / length(v) * vectorLength;
}

bool PointSet::isValid(const Site& site) const
{
	if (site.replicaIds.empty())
		return false;

	// Check if all replicas of the site have the same main replica
	const auto expectedMainReplica = getMainReplica(replicas[site.replicaIds[0]].vh->point());
	for (auto replicaId : site.replicaIds)
	{
		auto mainReplica = getMainReplica(replicas[replicaId].vh->point());
		if (length(expectedMainReplica - mainReplica) > 1e-10) // points aren't considered identical anymore
		{
			std::cout << "Found invalid site with the following main replicas of the points:\n";
			for (auto replicaId : site.replicaIds)
			{
				printCoordinates(getMainReplica(replicas[replicaId].vh->point()));
			}
			return false;
		}
	}
	
	return true;
}

PointSet::PointSet(int nPoints, int initType, double aspectRatio) : PointSet(nPoints, aspectRatio)
{
	switch (initType)
	{
	case 0:
		//initRandom(aspectRatio);
		fprintf(stderr, "Not yet implemented %d\n", initType);
		exit(1);
		break;
	case 1:
		//initDarts();
		fprintf(stderr, "Not yet implemented %d\n", initType);
		exit(1);
		break;
	case 2:
		//initJittered();
		fprintf(stderr, "Not yet implemented %d\n", initType);
		exit(1);
		break;
	case 3:
		//initGrid();
		fprintf(stderr, "Not yet implemented %d\n", initType);
		exit(1);
		break;
	default:
		fprintf(stderr, "Undefined initialization option %d\n", initType);
		exit(1);
	}
}

PointSet::PointSet(int nPoints, double* inputPoints, double aspectRatio) : PointSet(nPoints, aspectRatio)
{
	dts.resize(1);
	arrangements.resize(1);

	std::vector<ReplicaPrototype> replicaPrototypesTile;
	size_t dtId = 0;
	for (short replicaNumber = 0; replicaNumber < 9; ++replicaNumber) // tile1 replicated nine times in arrangement 0 to work on toroidal domain
	{
		replicaPrototypesTile.emplace_back(replicaNumber, dtId);
	}
	for (auto i = 0; i < n; ++i)
	{
		Point p(inputPoints[i], inputPoints[i + n]);
		addSite(p, replicaPrototypesTile);
	}
}

PointSet::PointSet(int nPoints, double* inputPoints, double* inputPoints2, double aspectRatio) : PointSet(nPoints, aspectRatio)
{
	twoTiles = true;
	dts.resize(3);
	arrangements.resize(3);

	std::vector<ReplicaPrototype> replicaPrototypesTile1;
	size_t dtId = 0;
	for (short replicaNumber = 0; replicaNumber < 9; ++replicaNumber) // tile1 replicated nine times in arrangement 0 to work on toroidal domain
	{
		replicaPrototypesTile1.emplace_back(replicaNumber, dtId);
	}
	dtId = 2;
	replicaPrototypesTile1.emplace_back(0, dtId); // tile1 in the center of arrangement 2 to enforce properties on tile1/tile2 border (0 is center tile)
	for (auto i = 0; i < n; ++i)
	{
		Point p(inputPoints[i], inputPoints[i + n]);
		addSite(p, replicaPrototypesTile1);
	}

	std::vector<ReplicaPrototype> replicaPrototypesTile2;
	dtId = 1;
	for (short replicaNumber = 0; replicaNumber < 9; ++replicaNumber) // tile2 replicated nine times in arrangement 1 to work on toroidal domain
	{
		replicaPrototypesTile2.emplace_back(replicaNumber, dtId);
	}
	dtId = 2;
	for (short replicaNumber = 1; replicaNumber < 9; ++replicaNumber) // tile2 in the eight outer tiles of arrangement 2 to enforce properties on tile1/tile2 border (1-8 are the outer tiles)
	{
		replicaPrototypesTile2.emplace_back(replicaNumber, dtId);
	}
	for (auto i = 0; i < n; ++i)
	{
		Point p(inputPoints2[i], inputPoints2[i + n]);
		addSite(p, replicaPrototypesTile2);
	}

	assert(replicas.size() == 3 * 9 * nPoints);
	assert(arrangements.size() == 3);
	assert(sites.size() == 2 * nPoints);
#ifndef NDEBUG
	// Check if size of data structures is as expected
	auto nSeventeenReplicas = 0; // sites that have nine replicas in the first and eight in the last arrangement
	auto nTenReplicas = 0; // sites that have nine replicas in the first and one in the last arrangement
	for (auto& site : sites)
	{
		auto nSiteReplicas = site.replicaIds.size();
		if (nSiteReplicas == 17)
			nSeventeenReplicas++;
		else
			if (nSiteReplicas == 10)
				nTenReplicas++;
	}
	assert(nSeventeenReplicas == nPoints);
	assert(nTenReplicas == nPoints);

	// Check if the points of all replicas have valid site ids
	for (auto& r : replicas)
	{
		auto siteId = r.vh->info().siteId;
		assert(0 <= siteId && siteId < sites.size());
	}

	// Check if all points in any of the dts have valid site ids
	for (auto& dt : dts)
	{
		for (auto vh = dt.finite_vertices_begin(); vh != dt.finite_vertices_end(); ++vh)
		{
			auto siteId = vh->info().siteId;
			assert(0 <= siteId && siteId < sites.size());
		}
	}

	// Check for all sites if the replicas have the same mainReplica
	for (auto i = 0; i < sites.size(); ++i)
	{
		//assert(isValid(sites[i]));
		isValid(sites[i]);
	}
	
#endif //NDEBUG
}

void PointSet::setDMin(double d)
{
	rel_dmin = d;
}

void PointSet::setRc(double r)
{
	rel_rc = r;
}

void PointSet::setSdA(double sd)
{
	sdA = sd;
}

int PointSet::ppo()
{
	auto iteration = 1; // todo: consider not counting the last iteration since it is only to check if all are stable

	allStable = false;
	for (auto& site : sites) site.isStable = false;

	for (; iteration <= maxIterations && !allStable; ++iteration)
	{
#ifdef PRINT_INFO
		std::cout << "Iteration" << iteration << std::endl;
#endif // PRINT_INFO
		for (auto& site : sites) site.becomeStable = true;
		allStable = true; // Assume the point set will be found stable. (This will be set to false as soon as a point is moved)
		for (const auto& arrangement : arrangements)
		{
			auto order = shuffle(arrangement.replicaIdsToIterate.size());
			for (size_t i = 0; i < arrangement.replicaIdsToIterate.size(); ++i)
			{
				const auto iRandom = order[i];
				const auto replicaId = arrangement.replicaIdsToIterate[iRandom];

				coverage(replicaId);
				conflict(replicaId);
				capacity(replicaId);
			}
		}

		for (auto& site : sites) site.isStable = site.becomeStable;

#ifdef PRINT_INFO
		std::cout << iteration << std::endl;
		minDistanceCheck();
#endif

#ifdef PLOT_POINTS
		plt::figure_size(1200, 780);
		std::vector<double> x;
		std::vector<double> y;
		getPoints(arrangements[0], x, y);
		plt::scatter(x, y);
		plt::xlim(0.0, ONE_X);
		plt::ylim(0.0, ONE_Y);
		plt::show();
#endif
	}

	return allStable ? iteration : -1;
}

void PointSet::minDistanceCheck()
{
	const auto dMin = rel_dmin * dHex;
	std::cout << "Allowed min edge length: " << dMin << "\n";
	std::cout << "The triangulations have the following edge lengths: \n";
	for (auto& dt : dts)
	{
		const auto mel = minEdgeLength(dt);
		std::cout << mel;
		if (mel < dMin)
			std::cout << " (dMin not satisfied)";
		std::cout << "\n";
	}
	std::cout.flush();
}

void PointSet::getPoints(double* outMatrix)
{
	if (twoTiles) // points of both point sets are returned in one matrix
	{
		auto arrangement = arrangements[0]; // contains the points of tile 0
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
		auto arrangement = arrangements[0]; // contains the points of tile 0
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

bool PointSet::isInMainReplica(const Point& p) const
{
	return 0 <= p.x() && p.x() < ONE_X && 0 <= p.y() && p.y() < ONE_Y;
}
