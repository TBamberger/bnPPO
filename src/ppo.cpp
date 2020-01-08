#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <complex>
#include <iostream>
#include <random>

#define FIXED_SEED // If this is defined a fixed seed is used for RNG to achive reproducible point sets; uses random seed otherwise

#define VL(x) sqrt((x).squared_length())

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                                  Point;
typedef K::Vector_2                                                 Vector;
typedef std::vector<Point>                                          Points;
typedef std::vector<Vector>                                         Vectors;
struct VInfo {                                                                  // Information stored in vertices in DT.
    int id;                                                                     // Index of point in t-map
};
struct FInfo {                                                                  // Information stored in faces in DT.
    Point c;                                                                    // Circumcenter, to avoid repeated calculation
};
typedef CGAL::Triangulation_vertex_base_with_info_2<VInfo, K>       Tvb;
typedef CGAL::Triangulation_face_base_with_info_2<FInfo, K>         Tfb;
typedef CGAL::Triangulation_data_structure_2<Tvb, Tfb>              Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      DT;
typedef DT::Vertex_handle                                           VH;
typedef DT::Face_handle                                             FH;
typedef DT::Vertex_circulator                                       VC;
typedef DT::Face_circulator                                         FC;

const double epsilon = 1e-12;

struct Statistics {                                                             // This structure is adapted from psa code "https://code.google.com/p/psa/"
    double mindist;                                                             // Global minimum nearest neighbor distance
    double avgmindist;                                                          // Average minimum nearest neighbor distance
    double orientorder;                                                         // Bond orientation order
    double effnyquist;                                                          // Effective Nyquist frequency
    double oscillations;                                                        // Oscillations, please refer to Heck et al. for meaning of these measures
    double coverageRadius;                                                      // We added this to monitor the size of holes
    double N4, N5, N6, N7, N8;                                                  // We added this to track number of neighbors percentage
    double sdA;                                                                 // Normalized standard deviation of areas of Voronoi cells
    double centroidDist;                                                        // Average distance of points to centroids of their cells
    double positionNorm;
    Statistics() : mindist(0), avgmindist(0), orientorder(0),
        effnyquist(0), oscillations(0), coverageRadius(0),
        N4(0), N5(0), N6(0), N7(0), N8(0) {};
};

inline bool isInRect(const Point& p, const Point& BL, const Point& TR) {
    return BL.x() <= p.x() && BL.y() <= p.y() &&
        p.x() < TR.x() && p.y() < TR.y();
}

inline double triangleType(Point& p1, Point& p2, Point& p3) {
    double sq1 = (p3 - p2).squared_length();
    double sq2 = (p1 - p3).squared_length();
    double sq3 = (p2 - p1).squared_length();
    if (sq1 < sq2) std::swap(sq1, sq2);
    if (sq1 < sq3) std::swap(sq1, sq3);
    return sq1 - (sq2 + sq3);                                                   // 0 for right, < 0 for acute, > 0 for obtuse
}

inline double triangleType(FC& fc) {
    return triangleType(
        fc->vertex(0)->point(),
        fc->vertex(1)->point(),
        fc->vertex(2)->point()
    );
}

class CPointSet
{
    std::default_random_engine re;                                              // Use this random engine for all rng related to the point set. This allows to have a fixed seeding per points set.
    std::uniform_real_distribution<> randX;                                     // Uniform random numbers over the whole domain in x direction
    std::uniform_real_distribution<> randY;                                     // Uniform random numbers over the whole domain in x direction
    int n;                                                                      // number of points in one period on an AABitmap
    double dhex;                                                                // Reference spacing of hexagonal packing.
    double rel_dmin = 0.87;                                                     // Relative minimum distance between points; twice the conflict radius
    double rel_rc = 0.65;                                                       // Relative maximum coverage radius.
    double sdA = 0.038600518;                                                   // Target standard deviation of cell areas; default is from Schlomer thesis p 64.
    bool allStable;                                                             // To implement termination criteria
    double ONE_X;                                                               // Make it possible to use another size for torroidal domain.
    double ONE_Y;
    double HALF_X;
    double HALF_Y;
    DT dt;                                                                      // Will maintain a Delaunay triangulation for relaxation and queries
    double maxShift;
    void updateFaceInfo() {
        DT::Finite_faces_iterator fit = dt.finite_faces_begin();
        for (; fit != dt.finite_faces_end(); fit++) {                           // Iterate through all (finite) faces in triangulation
            fit->info().c = dt.circumcenter(fit);                               // Circumcenter of face is one end of a Voronoi edge
        }
    }
    struct TSite {                                                              // This is to handle the Delaunay triangulation.
        Point p;
        VH vh[9];                                                               // Handles to up to 9 replicas to maintain toroidal domain
        Vector force;                                                           // Resultant attraction/repulsion force exerted on this point by neighbors
        bool isStable, becomeStable;                                            // Used during optimization to skip processing a point if neither it nor neighbors move.
        inline double x() { return p.x(); };
        inline double y() { return p.y(); };
    };
    std::vector<TSite> sites;                                                   // The final coordinates of points.
    double toroidallinearDistX(double x1, double x2) const {             // 1D Nearest distance between replicas of two points
        double dx = x1 - x2;                                                    // Find distance in primary period
        while (dx > HALF_X) dx -= ONE_X;                                            // If larger than half the period length another replica is certainly nearer
        while (dx < -HALF_X) dx += ONE_X;                                           // Same, but opposite ordering of points
        return dx;
    }
    double toroidallinearDistY(double x1, double x2) const {              // 1D Nearest distance between replicas of two points
        double dx = x1 - x2;                                                    // Find distance in primary period
        while (dx > HALF_Y) dx -= ONE_Y;                                            // If larger than half the period length another replica is certainly nearer
        while (dx < -HALF_Y) dx += ONE_Y;                                           // Same, but opposite ordering of points
        return dx;
    }
    double toroidalSqDist(Point& p1, Point& p2) const {                  // 2D Nearest distance; square to avoid costly square root operation
        double dx = toroidallinearDistX(p1.x(), p2.x());
        double dy = toroidallinearDistY(p1.y(), p2.y());
        return dx * dx + dy * dy;
    }
    double toroidalDist(Point& p1, Point& p2) const {
        return sqrt(toroidalSqDist(p1, p2));
    }
    Point mainReplica(Point& p) {
        double x = p.x(), y = p.y();
        while (x < 0) x += ONE_X;
        while (x >= ONE_X) x -= ONE_X;
        while (y < 0) y += ONE_Y;
        while (y >= ONE_Y) y -= ONE_Y;
        return Point(x, y);
    }

    Point replica(Point& p, int i) {                                     // Find one of the 9 replicas of a point
        i = (i + 4) % 9;                                                        // We make the middle replica at index 0
        double x = p.x() + (i % 3 - 1) * ONE_X;                                   // Add -ONE, 0, or ONE to x
        double y = p.y() + (i / 3 - 1) * ONE_Y;                                   // Same for y
        return Point(x, y);
    }
    Point marginBL, marginTR;                                                   // Points within these margins affect points in the main replica
    Point setSite(int index, Point p);
    void moveSite(int index, Point p);

    void moveSite(int index, Vector shift) {                             // shifts points relative to their current location
        if (shift.squared_length() > epsilon) {
            moveSite(index, sites[index].p + shift);
        }
    }
    // Relaxation forces:
    typedef Vector(CPointSet::* forceFunction) (int);                           // A pointer to a member function for calculating forces & return shift vectors
    Vector centroid(int index);
    Vector conflict(int index);
    Vector coverage(int index);
    Vector capacitySerial(int index);

    inline Point normalize(Point p) {                                           // Scale to unit & wrap around if points were shifted outside primary period
        p = mainReplica(p);
        double x = p.x() / ONE_X;
        double y = p.y() / ONE_Y;
        return Point(x, y);
    };
    void initRandom(double aspectRatio);
    void initDarts();
    void initJittered();
    void initGrid();
    unsigned* shuffle(const unsigned N);
    CPointSet(int number_of_points, double aspectRatio);
public:
    int stableCount;
    CPointSet(int number_of_points, int initType, double aspectRatio);                          // initType:   0: random, 1: darts, 2: jittered grid, 3: regular grid
    CPointSet(int nPoints, double* inputPoints, double aspectRatio); // todo: remove initType since it is not necessary here

    void getPoints(double* outMatrix); // outMatrix has to be allocated for 2*nPoints doubles (e.g. via mxCreateDoubleMatrix and mxGetPr)

    double PPO_serial(std::string seq); // seq specifies the order of local optimizations

    int getNumberOfPoints() { return n; }

    void setdmin(double d) { rel_dmin = d; };                                   // Set target NND for spring().
    void setRc(double r) { rel_rc = r; };
    void set_sdA(double sd) { sdA = sd; };
    bool isAllStable() { return allStable; };
    double getMaxShift() { return sqrt(maxShift); }
    Statistics GetStatistics();
    void setAllUnstable() {
        for (int i = 0; i < n; i++) sites[i].isStable = false;
        allStable = false;
    }
};

CPointSet::CPointSet(int number_of_points, double aspectRatio) :
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
    sites.resize(n);                                                          // The final location including shifts
}

CPointSet::CPointSet(int nPoints, int initType, double aspectRatio) : CPointSet(nPoints, aspectRatio)
{
    switch (initType)
    {
    case 0:
        initRandom(aspectRatio);
        break;
    case 1:
        initDarts();
        break;
    case 2:
        initJittered();
        break;
    case 3:
        initGrid();
        break;
    default:
        fprintf(stderr, "Undefined initialization option %d\n", initType);
        exit(1);
    }
}

CPointSet::CPointSet(int nPoints, double* inputPoints, double aspectRatio) : CPointSet(nPoints, aspectRatio)
{
    // todo: implement, process inputPoints
}

// @return: Randomly ordered list with the integers from 0 to n-1
unsigned* CPointSet::shuffle(const unsigned N)
{                
    auto randMax = std::uniform_int_distribution<>();
    unsigned* list = new unsigned[N];
    for (int i = 0; i < N; i++)
        list[i] = i;
    for (unsigned i = 0; i < N - 1; i++)
    {
        unsigned r = i + randMax(re) % (N - 1 - i);
        std::swap(list[i], list[r]);
    }
    return list;
}

void CPointSet::getPoints(double* outMatrix)
{
    for (int row = 0; row < n; row++)
    {
        Point p = mainReplica(sites[row].p);
        // If scaling point set back to [0,1)^2 is desired, use this line instead
        //Point p = normalize(sites[row].p); // Normalize and wrap back to unit torus

        outMatrix[row] = p.x();
        outMatrix[n + row] = p.y();
    }
}

void CPointSet::initRandom(double aspectRatio)
{
    for (auto i = 0; i < n; ++i)
    {
        Point p(randX(re), randY(re));
        setSite(i, p);
    }
}

/*
generate intial point set by Poisson-disk darts
*/
void CPointSet::initDarts()
{
    const int MAX = 10000; // To avoid infinite loop
    double r = 0.75 * dhex;
    double rr = r * r;
    for (int i = 0; i < n; i++) {
        Point p;
        bool accepted;
        do {
            p = Point(randX(re), randY(re));
            accepted = true;
            int counter = 0;
            for (int j = 0; (j < i) && accepted; j++) {
                accepted = accepted && (toroidalSqDist(p, sites[j].p) > rr);
            }
            if (counter++ >= MAX) {
                fprintf(stderr, "Darts failed at %d\n", i);
                exit(1);
            }
        } while (!accepted);
        setSite(i, p);
    }
}

/*
generate intial point set by jitter
*/
void CPointSet::initJittered()
{
    std::cerr << "jittered initialization not yet available for rectangular domains." << std::endl;
    exit(1);
    //int m = sqrt(n);
    //for (int i = 0; i < n; i++) {
    //    double x = i % m + 2 * drand48();
    //    double y = i / m + 2 * drand48();
    //    setSite(i, Point(x, y));
    //}
}

/*
generate intial point set regularly
*/
void CPointSet::initGrid()
{
    std::cerr << "grid initialization not yet available for rectangular domains." << std::endl;
    exit(1);
    //int m = sqrt(n);
    //for (int i = 0; i < n; i++) {
    //    double x = i % m + 0.5;
    //    double y = i / m + 0.5;
    //    setSite(i, Point(x, y));
    //}
}

Vector CPointSet::centroid(int index) {
    double a = 0, cx = 0, cy = 0;                                               // Cell area (actually twice the area) and centroid coordinates
    double XProduct;                                                            // Cross product of vectors to adjacent vertices
    FC fc = dt.incident_faces(sites[index].vh[0]), done(fc);
    do {
        Point p1 = dt.circumcenter(fc);
        Point p2 = dt.circumcenter(++fc);
        XProduct = p1.x() * p2.y() - p1.y() * p2.x();
        a += XProduct;                                                          // Accumulate areas
        cx += (p1.x() + p2.x()) * XProduct;
        cy += (p1.y() + p2.y()) * XProduct;
    } while (fc != done);
    cx /= 3.0 * a;
    cy /= 3.0 * a;
    return Point(cx, cy) - sites[index].p;                                      // Return shift from current position to centroid
};

Point CPointSet::setSite(int index, Point p) {                                  // Set location of the indexed point (in t-domain) and insert it in triangulation
    p = mainReplica(p);
    sites[index].p = p;                                                         // Save a handy copy of point coordinates
    sites[index].isStable = false;
    for (int i = 0; i < 9; i++) {                                               // We loop through the 9 replicas,
        //if (isInRect(replica(p, i), marginBL, marginTR)) {                      // if the location of a replica is within margin
        //    sites[index].vh[i] = dt.insert(replica(p, i));                      // insert replica in triangulation and keep handle to it
        //    sites[index].vh[i]->info().id = index;                              // Point the DT point back to map entry
        //}
        //else
        //{
        //    sites[index].vh[i] = NULL;
        //}

        // todo: is it justified to remove the margin to allow for other settings of ONE than sqrt(n)
        sites[index].vh[i] = dt.insert(replica(p, i));                      // insert replica in triangulation and keep handle to it
        sites[index].vh[i]->info().id = index;                              // Point the DT point back to map entry
    }
    return p;
};

void CPointSet::moveSite(int index, Point p) {                                  // Adjust location of indexed point (in t-domain) and update in triangulation
    double l = (p - sites[index].p).squared_length();
    maxShift = std::max(maxShift, l);
    sites[index].p = p;                                                         // Save a handy copy of updated point coordinates
    for (int i = 0; i < 9; i++) {
        if (sites[index].vh[i] != NULL)
            sites[index].vh[i] = dt.move(sites[index].vh[i], replica(p, i));
    }
    sites[index].becomeStable = false;                                          // Mark the point instable
    VC vc = dt.incident_vertices(sites[index].vh[0]), done(vc);
    do {                                                                        // Mark neighbors instable
        sites[vc->info().id].becomeStable = false;
    } while (++vc != done);
    allStable = false;                                                          // Mark the whole point set instable
};

Vector CPointSet::capacitySerial(int i)
{
    double d[20];    // Distance to neighbor (= 2 x distance to Voronoi edge)
    double el[20];   // length of the voronoi edges
    double area = 0; // Area of Voronoi cell
    FC fc2 = dt.incident_faces(sites[i].vh[0]), fc1(fc2++);                     // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
    VC vc = dt.incident_vertices(sites[i].vh[0], fc2), done(vc);                // The vertex sharing fc1 and fc2 with v[i].vh
    int m = 0;                                                                  // Number of neighbors
    Vector dir[20];                                                             // Direction vectors to neighbors
    int id[20];                                                                 // Id's of neighbors. We can't use the circulator for updating
    do {
        Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);             // Circumcenters of faces are endpoints of Voronoi cell edge
        area += c1.x() * c2.y() - c1.y() * c2.x();                              // Accumulate areas (rhs is cross product)
        el[m] = sqrt((c2 - c1).squared_length());                               // Length of Voronoi edge
        dir[m] = (vc->point() - sites[i].p);
        d[m] = sqrt(dir[m].squared_length());                                   
        dir[m] = dir[m] / d[m];                                                 // Normalize direction vector
        id[m] = vc->info().id;
        ++fc1;
        ++fc2;
        ++m;
    } while (++vc != done);
    area /= 2;
    double averageArea = ONE_X * ONE_Y / n;
    double dA = area - averageArea; // Required expansion or contraction
    if (fabs(dA) > sdA)
    {
        double sum_w = 0;
        for (int j = 0; j < m; j++)
            sum_w += el[j] * el[j];
        double pressure = -2 * dA / sum_w; // pressure per unit length of edges

        for (int j = 0; j < m; j++)
        {                                           
            Vector force = pressure * el[j] * dir[j];
            moveSite(id[j], force);
        }
    }
    return Vector(0, 0);
}

// test conflict resolution by pushing neighbors
Vector CPointSet::conflict(int index)
{
    double dmin = rel_dmin * dhex;
    Point& p = sites[index].p;
    bool conflict[30];
    Vector shift[30];
    int id[30];
    int m = 0;
    VC vc = dt.incident_vertices(sites[index].vh[0]), done(vc);
    do {
        Vector edge = vc->point() - p;
        double l = VL(edge);
        if (l < dmin) {
            conflict[m] = true;
            shift[m] = (1.001 * dmin / l - 1) * edge;
            id[m] = vc->info().id;
        }
        else conflict[m] = false;
        m++;
    } while (++vc != done);
    for (int i = 0; i < m; i++) {
        if (conflict[i]) moveSite(id[i], shift[i]);
    }
    return Vector(0, 0);
}

// A coverage routine which pulls neighbors
Vector CPointSet::coverage(int index)
{
    double rc = rel_rc * dhex;
    double dmin = rel_dmin * dhex;
    int id[30];
    double scale[30];
    Vector edge[30];
    int m = 0;
    FC fc = dt.incident_faces(sites[index].vh[0]), done(fc);
    VC vc = dt.incident_vertices(sites[index].vh[0], fc);
    do {
        vc++;
        edge[m] = sites[index].p - vc->point();
        id[m] = vc->info().id;
        if (triangleType(fc) <= 0) {
            Point c = dt.circumcenter(fc);
            double l = VL(c - sites[index].p);
            scale[m] = rc / l;
        }
        else scale[m] = 2; // > 1
        m++;
    } while (++fc != done);
    for (int i = 0; i < m; i++) {
        double scl = std::min(scale[i], scale[(i + 1) % m]);
        if (scl < 1) {
            Vector shift = (1 - scl) * edge[i];
            moveSite(id[i], shift);
        }
    }
    return Vector(0, 0);
}

double CPointSet::PPO_serial(std::string seq) {
    maxShift = 0;
    for (int i = 0; i < n; i++) sites[i].becomeStable = true;                   // Assume all points will be found stable.
    allStable = true;                                                           // Assume the point set will be found stable.
    unsigned* order = shuffle(n);
    for (int i = 0; i < n; i++) {                                               // Iterate through numbers up to maximum index
        int index = order[i];                                                   // Translate index according to the random order
        if (sites[index].isStable) continue;
        for (int k = 0; k < seq.length(); k++) {
            switch (seq[k]) {
            case '0': coverage(index); break;
            case '1': conflict(index); break;
            case '2': capacitySerial(index); break;
            }
        }
    }
    delete[] order;
    stableCount = 0;
    for (int i = 0; i < n; i++) {
        sites[i].isStable = sites[i].becomeStable;
        if (sites[i].becomeStable) stableCount++;
    }
    return sqrt(maxShift);                                                      // Return maximum shift which can be used to monitor convergence
}

Statistics CPointSet::GetStatistics() {                                         // Return statistics useful to monitor progress of optimization
    Statistics stats;                                                           // This function is adapted from PSA code.
    stats.mindist = std::max(ONE_X, ONE_Y);                                                        // The minimum distance can't go above ONE :)
    stats.avgmindist = 0;                                                       // Reset this to 0 to accumulate for averaging
    stats.orientorder = 0;                                                      // ~
    std::complex<double> acc = 0;
    unsigned long nacc = 0;
    for (int i = 0; i < n; i++) {
        double localmd = std::max(ONE_X, ONE_Y);
        std::complex<double> localacc = 0;
        VC vc = dt.incident_vertices(sites[i].vh[0]),
            done(vc), next;
        do {
            next = vc; ++next;
            const Point& v1 = vc->point();
            const Point& v2 = next->point();
            // Local mindist
            double dist = CGAL::squared_distance(sites[i].p, v1);
            localmd = std::min(localmd, (double)dist);
            // Orientational order
            std::complex<double> c1(v1.x(), v1.y());
            std::complex<double> c2(v2.x(), v2.y());
            localacc += std::polar(1.0, 6.0 * arg(c1 - c2));                    // This is how it's done in PSA (divide only at end); but I feel it needs review
            ++nacc;
        } while (++vc != done);
        stats.mindist = std::min(stats.mindist, localmd);
        stats.avgmindist += sqrtf(localmd);
        acc += abs(localacc);
    }
    // Coverage Radius
    double Rc = 0;
    Point BL(0, 0), TR(ONE_X, ONE_Y);                                               // We will consider only faces in one period
    DT::Finite_faces_iterator it;
    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++) {     // Iterate through all finite faces in triangulation
        if (isInRect(it->vertex(0)->point(), BL, TR)) {                         // It suffices to have one vertex within main period to consider the face
            Point c = dt.circumcenter(it);                                      // Farthest point in current face
            double r = CGAL::squared_distance(c, it->vertex(0)->point());       // Farthest distance in current face
            if (r > Rc) { Rc = r; }                                             // If larger than current farthest distance update the latter
        }
    }
    // Voronoi cell N-Gon types:                                   // Before orientation order this used to be a measure; see paper by Balzer et al.
    int histogram[10] = { 0,0,0,0,0,0,0,0,0,0 };                            // A simple histogram, where index 0 is interpreted as 3
    for (int i = 0; i < n; i++) {                                           // Iterate through all points
        VC vc = dt.incident_vertices(sites[i].vh[0]), done(vc);                 // We did not call neighbors to save the overhead of allocating Points
        int n = 0; do { n++; } while (++vc != done);                        // Count number of neighbors
        ++histogram[n - 3];                                                 // Increment the respective histogram bin; triangle is the smallest possible cell
    }
    // capacity variations
    updateFaceInfo();                                                       // Update circumcenters
    double aa = 0;
    for (int i = 0; i < n; i++) {                                           // First we need to calculate pressure of each cell
        double a = 0;                                                       // Area of Voronoi cell
        double XProduct;
        FC fc2 = dt.incident_faces(sites[i].vh[0]), fc1(fc2++), done(fc1);  // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
        do {
            Point& c1 = fc1->info().c, & c2 = fc2->info().c;                // Circumcenters of faces are endpoints of Voronoi cell edge
            XProduct = c1.x() * c2.y() - c1.y() * c2.x();
            a += XProduct;                                                  // Accumulate areas
            ++fc2;
        } while (++fc1 != done);
        a /= 2;
        aa += a * a;                                                        // We track convergence by monitoring variance in areas
    }
    double areaVariance = (aa / n) - 1;                                     // E(a^2) - (E(a)) ^ 2; where E(a) = 1;

    // position gradient; according to BNOT:
    double norm = 0;
    for (int i = 0; i < n; i++) {
        Vector dist = centroid(i);
        norm += dist.squared_length();
    }

    // Aggregate statistics:
    stats.mindist = sqrtf(stats.mindist) / dhex;                            // Nearest Neighbor Distance (NND) should
    stats.avgmindist /= n * dhex;                                           // be relative to dhex, the NND of hexagonal lattice
    stats.orientorder = abs(acc) / nacc;                                    // I don't quite get why we should average this once rather than average averages
    stats.coverageRadius = sqrtf(Rc) / dhex;                                // Size of largest hole (cf dmin); there should be an average too (cf davg)!
    stats.N4 = (double)histogram[1] / n;                                    // Ratios of various Voronoi cell types
    stats.N5 = (double)histogram[2] / n;                                    // 3 and 9+ are quite unlikely after relaxation
    stats.N6 = (double)histogram[3] / n;
    stats.N7 = (double)histogram[4] / n;
    stats.N8 = (double)histogram[5] / n;
    stats.sdA = sqrt(areaVariance);
    stats.positionNorm = sqrt(norm / n);
    return stats;
}

void optimizePattern(CPointSet ps, double *outMatrix, double aspectRatio)
{
	
}

void optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, int initType, double *outMatrix, double aspectRatio)
{
    // defaults from .bat file: dMin=0.85, rC=0.67, sda=0.02
    // defaults from code (have never used them): dMin=0.87, rC=0.65, sda=-1


    CPointSet ps(nPoints, initType, aspectRatio);
    ps.setdmin(dMin);
    ps.setRc(rC);
    double sdA = -1; // todo: areaDeltaMax is currently ignored
    if (sdA >= 0) ps.set_sdA(sdA);

    const auto iterations = 500; // pushPull bat file used 3500 iterations
    ps.setAllUnstable();
    for (int i = 0; i < iterations; ++i)
    {
        ps.PPO_serial("012");

        double stableRatio = (double)ps.stableCount / ps.getNumberOfPoints();
        fprintf(
            stderr,
            "%4d - stable = %6d, stable-ratio = %10.8f, max-step = %10.8f\n",
            i, ps.stableCount, stableRatio, ps.getMaxShift()
        );

        if (ps.isAllStable())
        {
            std::cout << "stable pattern reached at iteration " << i << std::endl;
            break;
        }
    }
	
    if (!ps.isAllStable())
        std::cout << "Did not reach a stable pattern in " << iterations << " iterations. Unfinished pattern is returned";
    ps.getPoints(outMatrix);
}
