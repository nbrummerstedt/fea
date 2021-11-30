#ifndef FEA_ELEMENT_QUADRILATERAL
#define FEA_ELEMENT_QUADRILATERAL

#include <array>
#include <vector>

#include "../euclid/Geometry"

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class Quad
{
	array<Point,4> nodes;

	public:
	Quad (	array<Point,4> p ) : nodes(p) {};

}; // class Quad

} // namespace fea

#endif
