#include "pointSet.h"

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

void PointSet::moveSite(size_t siteId, Point targetPoint)
{
	const auto currentPoint = getMainReplica(getPoint(sites[siteId].replicaIds[0])); // todo: getMainReplica könnte ich mir hier sparen wenn ich sicherstelle das id 0 immer in der mittleren Kachel ist
	targetPoint = getMainReplica(targetPoint);
	const Vector shift = targetPoint - currentPoint;

	for (auto replicaId : sites[siteId].replicaIds)
	{
		auto vh = replicas[replicaId].vh;
		auto newPosition = vh->point() + shift;
		replicas[replicaId].vh = dts[replicas[replicaId].dtId].move(vh, newPosition);

		// Mark neighbors unstable (has to be done for each replica since replicas might be in different triangulations
		// and therefore have different neighbors).
		VC vc = dts[replicas[replicaId].dtId].incident_vertices(vh), done(vc);
		do
		{
			if (!dts[replicas[replicaId].dtId].is_infinite(vc))
				sites[vc->info().id].becomeStable = false;
		}
		while (++vc != done);
	}
	sites[siteId].becomeStable = false;
	allStable = false;
}

void PointSet::moveSite(size_t siteId, const Vector shift)
{
	if (shift.squared_length() > epsilon)
	{
		moveSite(siteId, getPoint(sites[siteId].replicaIds[0]) + shift);
	}
}

void PointSet::coverage(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	const auto rc = rel_rc * dhex;
	int id[30];
	double scale[30];
	Vector edge[30];
	auto m = 0;

	FC fc = dt.incident_faces(replica.vh);
	const FC done(fc);
	VC vc = dt.incident_vertices(replica.vh, fc);
	do
	{
		++vc;
		edge[m] = p - vc->point();
		id[m] = vc->info().id;
		if (triangleType(fc) <= 0) // acute or right angled
		{
			Point c = dt.circumcenter(fc);
			double l = VL(c - p);
			scale[m] = rc / l;
		}
		else
		{
			scale[m] = 2; // > 1
		}
		m++;
	} while (++fc != done);

	for (auto i = 0; i < m; i++)
	{
		const double scl = std::min(scale[i], scale[(i + 1) % m]);
		if (scl < 1)
		{
			Vector shift = (1 - scl) * edge[i];
			moveSite(id[i], shift);
		}
	}
}

void PointSet::conflict(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	const double dMin = rel_dmin * dhex;
	bool conflict[30];
	Vector shift[30];
	int id[30];
	int m = 0;
	VC vc = dt.incident_vertices(replica.vh);
	const VC done(vc);
	do
	{
		Vector edge = vc->point() - p;
		double l = VL(edge);
		if (l < dMin) {
			conflict[m] = true;
			shift[m] = (1.001 * dMin / l - 1) * edge;
			id[m] = vc->info().id;
		}
		else conflict[m] = false;
		m++;
	} while (++vc != done);

	for (auto i = 0; i < m; i++)
	{
		if (conflict[i]) moveSite(id[i], shift[i]);
	}
}

void PointSet::capacity(size_t replicaId)
{
	auto& replica = replicas[replicaId];
	auto& dt = dts[replica.dtId];
	auto& p = getPoint(replicaId);

	double d[20];    // Distance to neighbor (= 2 x distance to Voronoi edge)
	double el[20];   // length of the voronoi edges
	double area = 0; // Area of Voronoi cell
	FC fc2 = dt.incident_faces(replica.vh), fc1(fc2++);               // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
	VC vc = dt.incident_vertices(replica.vh, fc2);                    // The vertex sharing fc1 and fc2 with v[i].vh
	const VC done(vc);
	int m = 0;                                                        // Number of neighbors
	Vector dir[20];                                                   // Direction vectors to neighbors
	int id[20];                                                       // Id's of neighbors. We can't use the circulator for updating
	do
	{
		Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);   // Circumcenters of faces are endpoints of Voronoi cell edge
		area += c1.x() * c2.y() - c1.y() * c2.x();                    // Accumulate areas (rhs is cross product)
		el[m] = sqrt((c2 - c1).squared_length());                     // Length of Voronoi edge
		dir[m] = (vc->point() - p);
		d[m] = sqrt(dir[m].squared_length());
		dir[m] = dir[m] / d[m];                                       // Normalize direction vector
		id[m] = vc->info().id;
		++fc1;
		++fc2;
		++m;
	}
	while (++vc != done);
	area /= 2;
	const double averageArea = ONE_X * ONE_Y / n;
	const double dA = area - averageArea; // Required expansion or contraction
	if (fabs(dA) > sdA)
	{
		double sum_w = 0;
		for (int j = 0; j < m; j++)
		{
			sum_w += el[j] * el[j];
		}
		const double pressure = -2 * dA / sum_w; // pressure per unit length of edges
		for (auto j = 0; j < m; j++)
		{
			const Vector force = pressure * el[j] * dir[j];
			moveSite(id[j], force);
		}
	}
}

PointSet::PointSet(int number_of_points, double aspectRatio) :
#ifdef FIXED_SEED
	re(12345)
#else
	re(std::random_device{}())
#endif // FIXED_SEED
{
	n = number_of_points;                                                     // Number of points in one period
	ONE_X = aspectRatio;
	ONE_Y = 1;
	HALF_X = 0.5 * ONE_X;
	HALF_Y = 0.5 * ONE_Y;
	randX = std::uniform_real_distribution<>(0, ONE_X);
	randY = std::uniform_real_distribution<>(0, ONE_Y);
	dhex = sqrt(ONE_X * ONE_Y * 2 / (sqrt(3) * n));                           // Maximum packing distance
	double margin = std::min(10 / sqrt(n), 1.0);                              // Margin for toroidal domain. Heuristically, use 10 layers of points.
	marginBL = Point(-margin * ONE_X, -margin * ONE_Y);                       // Bottom-left of primary period + margin
	marginTR = Point((1 + margin) * ONE_X, (1 + margin) * ONE_Y);             // Top-right. In our convention BL is included, TR is excluded
	//sites.resize(2 * n);                                                          // The final location including shifts
}

PointSet::PointSet(int nPoints, int initType, double aspectRatio) : PointSet(nPoints, aspectRatio)
{
	switch (initType)
	{
	case 0:
		//initRandom(aspectRatio);
		fprintf(stderr, "Undefined initialization option %d\n", initType);
		exit(1);
		break;
	case 1:
		//initDarts();
		fprintf(stderr, "Undefined initialization option %d\n", initType);
		exit(1);
		break;
	case 2:
		//initJittered();
		fprintf(stderr, "Undefined initialization option %d\n", initType);
		exit(1);
		break;
	case 3:
		//initGrid();
		fprintf(stderr, "Undefined initialization option %d\n", initType);
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
	
	for (auto i = 0; i < n; ++i)
	{
		Point p(inputPoints[i], inputPoints[i + n]);

		Site site;
		const auto siteId = sites.size();
		for (auto j = 0; j < 9; ++j)
		{
			Replica r;
			const auto replicaId = replicas.size();
			r.vh = dts[0].insert(createReplica(p, j));
			r.vh->info().id = siteId;
			r.dtId = 0;
			replicas.push_back(r);
			site.replicaIds.push_back(replicaId);
			if (j == 0) // replicas in center tile
			{
				arrangements[0].replicaIdsToIterate.push_back(replicaId);
			}
		}
		sites.push_back(site);
	}
}

PointSet::PointSet(int nPoints, double* inputPoints, double* inputPoints2, double aspectRatio) : PointSet(nPoints, aspectRatio) // todo: 
{
	//dts.resize(3);

	//// toroidal tile 1
	//for (auto i = 0; i < n; ++i)
	//{
	//	Point p(inputPoints[i], inputPoints[i + n]);

	//	Site site;
	//	const auto siteId = sites.size();
	//	for (auto j = 0; j < 9; ++j)
	//	{
	//		Replica r;
	//		const auto replicaId = replicas.size();

	//		r.vh = dts[0].insert(createReplica(p, j));
	//		r.vh->info().id = siteId;
	//		r.dtId = 0;
	//		replicas.push_back(r);
	//		site.replicaIds.push_back(replicaId);
	//	}
	//	sites.push_back(site);
	//}

	//// toroidal tile 2
	//for (auto i = 0; i < n; ++i)
	//{
	//	Point p(inputPoints2[i], inputPoints2[i + n]);

	//	Site site;
	//	const auto siteId = sites.size();
	//	for (auto j = 0; j < 9; ++j)
	//	{
	//		Replica r;
	//		const auto replicaId = replicas.size();

	//		r.vh = dts[1].insert(createReplica(p, j));
	//		r.vh->info().id = siteId;
	//		r.dtId = 1;
	//		replicas.push_back(r);
	//		site.replicaIds.push_back(replicaId);
	//	}
	//	sites.push_back(site);
	//}

	//// tile1/tile2 border
	//for (auto i = 0; i < n; ++i) // place tiles of tile 1 in the center tile
	//{
	//	Point p(inputPoints[i], inputPoints[i + n]);

	//	Site site;
	//	const auto siteId = sites.size();
	//	for (auto j = 0; j < 9; ++j)
	//	{
	//		Replica r;
	//		const auto replicaId = replicas.size();

	//		r.vh = dts[3].insert(createReplica(p, j));
	//		r.vh->info().id = siteId;
	//		r.dtId = 3;
	//		replicas.push_back(r);
	//		site.replicaIds.push_back(replicaId);
	//	}
	//	sites.push_back(site);
	//}

	//for (auto i = 0; i < n; ++i)
	//{
	//	Point p(inputPoints2[i], inputPoints2[i + n]);

	//	Site site;
	//	const auto siteId = sites.size();
	//	for (auto j = 0; j < 9; ++j)
	//	{
	//		Replica r;
	//		const auto replicaId = replicas.size();

	//		r.vh = dts[3].insert(createReplica(p, j));
	//		r.vh->info().id = siteId;
	//		r.dtId = 3;
	//		replicas.push_back(r);
	//		site.replicaIds.push_back(replicaId);
	//	}
	//	sites.push_back(site);
	//}
}

int PointSet::ppo()
{
	allStable = false;
	for (auto& site : sites) site.isStable = false;

	auto iteration = 1; // todo: consider not counting the last iteration since it is only to check if all are stable
	for (; iteration <= maxIterations && !allStable; ++iteration)
	{
		for (auto& site : sites) site.becomeStable = true;
		allStable = true; // Assume the point set will be found stable. (This will be set to false as soon as a point is moved)
		for (const auto& arrangement : arrangements)
		{
			for (auto replicaId : arrangement.replicaIdsToIterate) // todo: consider traversing the points in site in random order
			{
				coverage(replicaId);
				conflict(replicaId);
				capacity(replicaId);
			}
		}
		for (auto& site : sites) site.isStable = site.becomeStable;
		std::cout << iteration << std::endl;
	}
	return allStable ? iteration : -1;
}
