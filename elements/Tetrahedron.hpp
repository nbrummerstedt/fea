#ifndef FEA_ELEMENT_TETRAHEDRON
#define FEA_ELEMENT_TETRAHEDRON

#include <array>
#include <vector>

namespace fea {

using std::array;
using std::vector;
using Euclid::Point;

class Tet
{
	array<Point,4>  global_coordinates;

	public:
	Tet ( array<Point,4> p ) : global_coordinates(p)	{};
	Volume( const array<Point,4> & );
	Shape ( const Point & );
	Shape ( const double &, const double &, const double & );
	Shape ( const IntegrationPoint & );

	private :
	// integration_points_tetrahedra_barycentric
	// From Zienkiewicz CH6, using volume coordinates
	// Conversion to barycentric (volume) coordinates from intrinsic rst-coordinates
	// L1 = 1-r-s-t, L2 = r, L3 = s, L4 = t
	// See Cook et al. p 266, and (less readable) Zienkiewicz section 6.3
	static constexpr double
	integration_points_tetrahedra_barycentric[5][5]
	{//   L1   L2   L3   L4   W
		{ 1/4, 1/4, 1/4, 1/4,-4/5 }, // ip1
		{ 1/2, 1/6, 1/6, 1/6, 9/20}, // ip2
		{ 1/6, 1/2, 1/6, 1/6, 9/20}, // ip3
		{ 1/6, 1/6, 1/2, 1/6, 9/20}, // ip4
		{ 1/6, 1/6, 1/6, 1/2, 9/20}  // ip5
	};

	// integration_points_intrinsic_4
	// Integration points in intrinsic (natural) coordinates
	// as used by Cook et al. Table 7.4-2
	// Can also be found here:
	// https://www.cfd-online.com/Wiki/Code:_Quadrature_on_Tetrahedra
	static constexpr array<FEA::IntegrationPoint,4>
	integration_points_intrinsic_4 {{
		{ 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.25/6 },
		{ 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.25/6 },
		{ 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.25/6 },
		{ 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.25/6 }
	}};

	// integration_points_intrinsic_5
	// Common integration scheme with negative weight center point
	private : static constexpr array<FEA::IntegrationPoint,5>
	integration_points_intrinsic_5 {{
		{ 1/4., 1/4., 1/4.,-4/5. },
		{ 1/6., 1/6., 1/6., 9/20.},
		{ 1/2., 1/6., 1/6., 9/20.},
		{ 1/6., 1/2., 1/6., 9/20.},
		{ 1/6., 1/6., 1/2., 9/20.}
	}};

}; // class Tet


// Volume
// Computes the volume of a Tet
static constexpr double
Tet::Volume( std::array<Euclid::Point,4> T ) {
	double Vp { (T[3].x()-T[0].x()) * (
					  (T[1].y()-T[0].y())*(T[2].z()-T[0].z())
					- (T[1].z()-T[0].z())*(T[2].y()-T[0].y())
				)
				+ (T[3].y()-T[0].y()) * (
					  (T[1].z()-T[0].z())*(T[2].x()-T[0].x())
					- (T[1].x()-T[0].x())*(T[2].z()-T[0].z())
				)
				+ (T[3].z()-T[0].z()) * (
					  (T[1].x()-T[0].x())*(T[2].y()-T[0].y())
					- (T[1].y()-T[0].y())*(T[2].x()-T[0].x())
				) };
	return Vp/6.;
};


// Shape
// Returns N = 4 shape functions (scalar) for a linear tetrahedron
// evaluated at local coordinate (r,s,t)
constexpr array<double,4>
Tet::Shape ( Euclid::Point & p )
{
	return array<double,4> {{1-p.x()-p.y()-p.z(),p.x(),p.y(),p.z()}};
};
constexpr array<double,4>
Tet::Shape ( double & r, double & s, double & t )
{
	return array<double,4> {{1-r-s-t,r,s,t}};
};
constexpr array<double,4>
Tet::Shape ( const FEA::IntegrationPoint & p )
{
	return array<double,4> {{1-p.x-p.y-p.z,p.x,p.y,p.z}};
};

// IntegrationPoints
// Returns a vector of integration points,
// each given in an standard-array with values {xi,eta,zeta,weight}
const vector<FEA::IntegrationPoint>
Tet::IntegrationPoints ( unsigned int count = 4 )
{
	vector<FEA::IntegrationPoint> ip;
	if (count==4) { for ( auto p : integration_points_intrinsic_4 ) ip.push_back(p); return ip; };
	if (count==5) { for ( auto p : integration_points_intrinsic_4 ) ip.push_back(p); return ip; };
	throw "Invalid input!";
};

// levelset_to_cutcase
// Computes the cut case given a level value for the isoline to cut with
// Usually the cut will be performed at a cut level of zero.
constexpr int
Tet::levelset_to_cutcase ( array<double,4> level_set , double cut_level )
{
	int CaseIndex { 0 };
	for (int i = 0; i < 4; i++)
		CaseIndex |= ( level_set[i] < cut_level ) << i;
	return CaseIndex;
};

} // namespace fea

#endif
