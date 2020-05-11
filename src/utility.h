#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <vector>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using Vector = K::Vector_2;
using Points = std::vector<Point>;
using Vectors = std::vector<Vector>;
struct VertexInfo
{
	int siteId;
};

using Tvb = CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, K>;
using Tfb = CGAL::Triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Tvb, Tfb>;
using DT = CGAL::Delaunay_triangulation_2<K, Tds>;
using VH = DT::Vertex_handle;
using FH = DT::Face_handle;
using VC = DT::Vertex_circulator;
using FC = DT::Face_circulator;

const double epsilon = 1e-12;

inline double length(const Vector& v)
{
	return sqrt(v.squared_length());
}

inline bool isInRect(const Point& p, const Point& bl, const Point& tr) // todo: test
{
	return bl.x() <= p.x() && bl.y() <= p.y() &&
		p.x() < tr.x() && p.y() < tr.y();
}

inline double triangleType(Point& p1, Point& p2, Point& p3) // todo: change return type to enum
{
	auto sq1 = (p3 - p2).squared_length();
	auto sq2 = (p1 - p3).squared_length();
	auto sq3 = (p2 - p1).squared_length();
	if (sq1 < sq2) std::swap(sq1, sq2);
	if (sq1 < sq3) std::swap(sq1, sq3);
	return sq1 - (sq2 + sq3); // 0 for right, < 0 for acute, > 0 for obtuse
}

inline double triangleType(FC& fc) // todo: return enum
{
	return triangleType(
		fc->vertex(0)->point(),
		fc->vertex(1)->point(),
		fc->vertex(2)->point()
	);
}

/// @return: Minimum edge length in the dt
inline double minEdgeLength(DT& dt)
{
	auto minEdgeLength = std::numeric_limits<double>::max();
	for (const auto& edge : dt.finite_edges())
	{
		Point p1 = edge.first->vertex(DT::cw(edge.second))->point();
		Point p2 = edge.first->vertex(DT::ccw(edge.second))->point();
		const double edgeLength = std::sqrt((p1 - p2).squared_length());
		if (edgeLength < minEdgeLength)
			minEdgeLength = edgeLength;
	}
	return minEdgeLength;
}

inline void printCoordinates(Point p)
{
	std::cout << p.x() << ", " << p.y() << std::endl;
}
