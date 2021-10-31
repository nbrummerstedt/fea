#ifndef FEA_CUTFEM_HEXAHEDRON
#define FEA_CUTFEM_HEXAHEDRON

#include "../elements/Hexahedron.hpp"

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class HexCut
{
	array<double,8> levels;
	public:
	HexCut ( array<Point,8> p , array<double,8> v ) : Hex(p), levels(v) {};

}; // class HexCut

} // namespace FEA

#endif
