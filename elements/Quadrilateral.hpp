#ifndef FEA_ELEMENT_QUADRILATERAL
#define FEA_ELEMENT_QUADRILATERAL

#include <array>
#include <vector>

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class Quad
{
	array<Point,4> 	global_coordinates;

	public:
	Quad (	array<Point,4> p ) : global_coordinates(p) {};

}; // end class Quad

} // end namespace fea

#endif
