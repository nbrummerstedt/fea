#ifndef FEA_ELEMENT_TRIANGLE
#define FEA_ELEMENT_TRIANGLE

#include <array>
#include <vector>

#include "../euclid/Geometry"

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class Tri
{
	array<Point,3> nodes;

	public:
	Tri ( array<Point,3> p ) : nodes(p)	{};

	private:
	// integration_points_triangle_intrinsic
	// Intrinsic (natural) coordinate formulation
	static constexpr double integration_points_intrinsic[4][3]
	{ //  r    s    w
		{ 1/3, 1/3,-27/96 },
		{ 1/5, 1/5, 25/96 },
		{ 3/5, 1/5, 25/96 },
		{ 1/5, 3/5, 25/96 }
	};

}; // class Tri

} // namespace fea


#endif
