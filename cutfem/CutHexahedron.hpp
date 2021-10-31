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
	static constexpr int 			levelset_to_cutcase 	( array<double,8> level_set , double cut_level );
	static constexpr int 			cutcase_to_basecase 	( int cutcase );
	static constexpr int 			cutcase_to_permutation 	( int cutcase );
	static constexpr array<bool,12> basecase_to_edges 		( int basecase );
	static constexpr Subdivision 	basecase_to_subdivision ( int basecase );

	// CutFEM related methods that return vectors
	const vector<Point> 			CutCoordinatesLocal        	( int d_i = -1, double d_l = 0. );
	const vector<Point> 			SubVertexCoordinatesLocal  	( int d_i = -1, double d_l = 0. );
	const vector<Point>  			SubVertexCoordinatesGlobal 	( int d_i = -1, double d_l = 0. );
	const vector<bool> 				SubElementPhases 			( );
	const vector<array<int,4>> 		SubElementIndices 			( );
	const vector<array<Point,4>> 	SubElementVertices 			( int d_i = -1, double d_l = 0. );

}; // class HexCut



} // namespace FEA

#endif
