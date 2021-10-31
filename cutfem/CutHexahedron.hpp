#ifndef FEA_CUTFEM_HEXAHEDRON
#define FEA_CUTFEM_HEXAHEDRON

#include "../elements/Hexahedron.hpp"

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class HexCut : public Hex
{
	array<double,8> levels;
	public:
	HexCut ( array<Point,8> p , array<double,8> v ) : Hex(p), levels(v) {};

	const vector<IntegrationPoint> SubIntegrationPoints	( int d_i = -1, double d_l = 0. );

	private:

	// Internal type
	typedef struct {
		const unsigned vertex_count;
		const unsigned element_count;
		const array<array<int,4>,22> connectivity;
	} Subdivision;

	// static data expressed as constexpr methods
	static constexpr int 			cutcase 	( array<double,8> level_set , double cut_level );
	static constexpr int 			basecase 	( int cutcase );
	static constexpr int 			permutation 	( int cutcase );
	static constexpr array<bool,12> edges 		( int basecase );
	static constexpr Subdivision 	subdivide ( int basecase );

	// CutFEM related methods that return vectors
	const vector<Point> 			CutCoordinatesLocal        	( int d_i = -1, double d_l = 0. );
	const vector<Point> 			SubVertexCoordinatesLocal  	( int d_i = -1, double d_l = 0. );
	const vector<Point>  			SubVertexCoordinatesGlobal 	( int d_i = -1, double d_l = 0. );
	const vector<bool> 				SubElementPhases 			( );
	const vector<array<int,4>> 		SubElementIndices 			( );
	const vector<array<Point,4>> 	SubElementVertices 			( int d_i = -1, double d_l = 0. );

	static const array<int,256> cutcase_to_basecase {{
		0,  1,  1,  2,  1,  3,  2,  5,  1,  2,  3,  5,  2,  5,  5,  8,
		1,  2,  3,  5,  4,  6,  6, 11,  3,  5,  7,  9,  6, 11, 12,  5,
		1,  3,  2,  5,  3,  7,  5,  9,  4,  6,  6, 11,  6, 12, 11,  5,
		2,  5,  5,  8,  6, 12, 11,  5,  6, 11, 12,  5, 10,  6,  6,  2,
		1,  4,  3,  6,  2,  6,  5, 11,  3,  6,  7, 12,  5, 11,  9,  5,
		3,  6,  7, 12,  6, 10, 12,  6,  7, 12, 13,  7, 12,  6,  7,  3,
		2,  6,  5, 11,  5, 12,  8,  5,  6, 10, 12,  6, 11,  6,  5,  2,
		5, 11,  9,  5, 11,  6,  5,  2, 12,  6,  7,  3,  6,  4,  3,  1,
		1,  3,  4,  6,  3,  7,  6, 12,  2,  5,  6, 11,  5,  9, 11,  5,
		2,  5,  6, 11,  6, 12, 10,  6,  5,  8, 12,  5, 11,  5,  6,  2,
		3,  7,  6, 12,  7, 13, 12,  7,  6, 12, 10,  6, 12,  7,  6,  3,
		5,  9, 11,  5, 12,  7,  6,  3, 11,  5,  6,  2,  6,  3,  4,  1,
		2,  6,  6, 10,  5, 12, 11,  6,  5, 11, 12,  6,  8,  5,  5,  2,
		5, 11, 12,  6, 11,  6,  6,  4,  9,  5,  7,  3,  5,  2,  3,  1,
		5, 12, 11,  6,  9,  7,  5,  3, 11,  6,  6,  4,  5,  3,  2,  1,
		8,  5,  5,  2,  5,  3,  2,  1,  5,  2,  3,  1,  2,  1,  1,  0
	}};

}; // class HexCut



} // namespace FEA

#endif
