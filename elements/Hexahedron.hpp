#ifndef ARGYRIS_ELEMENT_HEXAHEDRON
#define ARGYRIS_ELEMENT_HEXAHEDRON

#include <vector>
#include <array>
#include <iostream>
#include <fstream>

namespace Argyris {

using std::array;
using std::vector;
using Euclid::Point;

class Hexahedron
{
	protected: array<Point,8> coordinates;
	public:
		Hex ( array<double,8> x, array<double,8> y, array<double,8> z ) : coordinates(doubles_to_points(x,y,z)) {};
		const vector<IntegrationPoint> 	IntegrationPoints 	( unsigned count = 8 );

	protected:

		// Coordinate extraction
		constexpr array<double,8> x () const;
		constexpr array<double,8> y () const;
		constexpr array<double,8> z () const;

		// FEM Standard functions
		static const array<double,8> 			Shape 			( Point & );
		static const Eigen::Matrix<double,3,8> 	ShapeDerivative ( Point & );
		static const Eigen::Matrix<double,3,3> 	JacobianMatrix  ( Point & , array<Point,8> & );
		static const double  					Jacobian 		( Point & , array<Point,8> & );
		static const double  					Jacobian 		( IntegrationPoint & , array<Point,8> & );

		// This-referring FEM methods
		const Eigen::Matrix<double,6,24>   		StrainDisplacementMatrix ( IntegrationPoint & , array<Point,8> & );

		// Useful lookups
		static constexpr double					Volume					( array<Point,8> & );
		const int 								PermuteNode				( int , int );
		const array<Point,8> 					PermuteCoordinates 		( array<Point,8> , int );

		// Almost data tables
		static constexpr array<Point,8> 		doubles_to_points 		( array<double,8>, array<double,8>, array<double,8> );
		static constexpr array<int,2> 			edge_to_endpoints 		( int );
		static constexpr array<int,8> 			permutation_to_indices 	( int );

		// FEM data tables
		// const (not constexpr) simply so they don't need initialisers here
		static const array<Point,8> 			local_coordinates;
		static const array<IntegrationPoint,8>  integration_points_linear;
		static const array<IntegrationPoint,27> integration_points_quadratic;
		static const array<array<int,8>,48> 	permutation_list;

}; // class Hex

// IntegrationPoints
// Returns a vector of standard integration array<Point,8> for a Hex,
// each given in an standard-array with values {xi,eta,zeta,weight}
inline const vector<IntegrationPoint>
Hex::IntegrationPoints ( unsigned int count )
{
	vector<IntegrationPoint> ips;
	if (count==8 ) { for ( auto p : integration_points_linear    ) ips.push_back(p); return ips; };
	if (count==27) { for ( auto p : integration_points_quadratic ) ips.push_back(p); return ips; };
	throw "Invalid input!";
};

// Coordinate extraction
inline constexpr array<double,8>
Hex::x () const { return array<double,8>{{
	coordinates[0].x(), coordinates[1].x(),
	coordinates[2].x(), coordinates[3].x(),
	coordinates[4].x(), coordinates[5].x(),
	coordinates[6].x(), coordinates[7].x() }};
};
inline constexpr array<double,8>
Hex::y () const { return array<double,8> {
	coordinates[0].y(), coordinates[1].y(),
	coordinates[2].y(), coordinates[3].y(),
	coordinates[4].y(), coordinates[5].y(),
	coordinates[6].y(), coordinates[7].y() };
};
inline constexpr array<double,8>
Hex::z () const { return array<double,8> {
	coordinates[0].z(), coordinates[1].z(),
	coordinates[2].z(), coordinates[3].z(),
	coordinates[4].z(), coordinates[5].z(),
	coordinates[6].z(), coordinates[7].z() };
};


// Shape
inline const array<double,8>
Hex::Shape ( Point & v )
{
	double x {v.x()};
	double y {v.y()};
	double z {v.z()};
	return {{ 	0.125*(1.-x)*(1.-y)*(1.-z),
				0.125*(1.+x)*(1.-y)*(1.-z),
				0.125*(1.+x)*(1.+y)*(1.-z),
				0.125*(1.-x)*(1.+y)*(1.-z),
				0.125*(1.-x)*(1.-y)*(1.+z),
				0.125*(1.+x)*(1.-y)*(1.+z),
				0.125*(1.+x)*(1.+y)*(1.+z),
				0.125*(1.-x)*(1.+y)*(1.+z)  }};
};

// shapeDerivative
// Row N (1-3) contains all eight shape functions
// derived with respect to direction N
inline const
Eigen::Matrix<double,3,8>
Hex::ShapeDerivative ( Point & v ) {
	double x {v.x()}; double y {v.y()}; double z {v.z()};
	Eigen::Matrix<double,3,8> D;
	D << -(1.-y)*(1.-z), (1.-y)*(1.-z), (1.+y)*(1.-z),-(1.+y)*(1.-z),-(1.-y)*(1.+z), (1.-y)*(1.+z), (1.+y)*(1.+z),-(1.+y)*(1.+z),
		 -(1.-x)*(1.-z),-(1.+x)*(1.-z), (1.+x)*(1.-z), (1.-x)*(1.-z),-(1.-x)*(1.+z),-(1.+x)*(1.+z), (1.+x)*(1.+z), (1.-x)*(1.+z),
		 -(1.-x)*(1.-y),-(1.+x)*(1.-y),-(1.+x)*(1.+y),-(1.-x)*(1.+y), (1.-x)*(1.-y), (1.+x)*(1.-y), (1.+x)*(1.+y), (1.-x)*(1.+y);
	return D * 0.125;
};

// jacobianMatrix
inline const Eigen::Matrix<double,3,3>
Hex::JacobianMatrix ( Point & vl , array<Point,8> & ng )
{
	Eigen::Matrix<double,8,3> XG;
	XG << 	ng[0].x(), ng[0].y(), ng[0].z(),
			ng[1].x(), ng[1].y(), ng[1].z(),
			ng[2].x(), ng[2].y(), ng[2].z(),
			ng[3].x(), ng[3].y(), ng[3].z(),
			ng[4].x(), ng[4].y(), ng[4].z(),
			ng[5].x(), ng[5].y(), ng[5].z(),
			ng[6].x(), ng[6].y(), ng[6].z(),
			ng[7].x(), ng[7].y(), ng[7].z();
	Eigen::Matrix<double,3,8> D { Hex::ShapeDerivative( vl ) };
	Eigen::Matrix<double,3,3> J = D * XG;
	return J;
};

// Jacobian
// Determinant to jacobi-matrix
inline const double
Hex::Jacobian ( Point & p , array<Point,8> & nodes_global )  {
	Eigen::Matrix<double,3,3> J { JacobianMatrix ( p , nodes_global ) };
	return { J.determinant() };
};
inline const double
Hex::Jacobian ( IntegrationPoint & ip , array<Point,8> & nodes_global )  {
	Point p {ip.x,ip.y,ip.z};
	Eigen::Matrix<double,3,3> J { JacobianMatrix ( p , nodes_global ) };
	return { J.determinant() };
};

// Strain-Displacement matrix
inline const Eigen::Matrix<double,6,24>
Hex::StrainDisplacementMatrix ( IntegrationPoint & evaluate_at , array<Point,8> & nodes_global )
{
	// Inverse Jacobian Matrix at evaluation point expanded to three dimensions
	Point p {evaluate_at.x,evaluate_at.y,evaluate_at.z};
	Eigen::Matrix<double,3,3> J { JacobianMatrix(p,nodes_global) };
	Eigen::Matrix<double,3,3> g { J.inverse() };
	Eigen::Matrix<double,9,9> G = Eigen::Matrix<double,9,9>::Zero();
	for (int i = 0; i < 3*3; i=i+3 ) { G.block<3,3>(i,i) = g; }
	Eigen::Matrix<double,6,9> L;
	L << 	1,0,0, 0,0,0, 0,0,0,
			0,0,0, 0,1,0, 0,0,0,
			0,0,0, 0,0,0, 0,0,1,
			0,1,0, 1,0,0, 0,0,0,
			0,0,0, 0,0,1, 0,1,0,
			0,0,1, 0,0,0, 1,0,0;
	Eigen::Matrix<double,3,8> d { ShapeDerivative ( p ) };
	Eigen::Matrix<double,9,24> D = Eigen::Matrix<double,9,24>::Zero();
	int s;
	for (int i = 1; i <= 8; i++ ) {
		s = (i*3)-2-1;
		D(0,s+0) = d(0,i-1);
		D(1,s+0) = d(1,i-1);
		D(2,s+0) = d(2,i-1);
		D(3,s+1) = d(0,i-1);
		D(4,s+1) = d(1,i-1);
		D(5,s+1) = d(2,i-1);
		D(6,s+2) = d(0,i-1);
		D(7,s+2) = d(1,i-1);
		D(8,s+2) = d(2,i-1);
	}
	Eigen::Matrix<double,6,24> B = L*G*D;
	return B;
}

// Volume
// Computes the volume of any Hex
inline constexpr double
Hex::Volume ( array<Point,8> & vertices )
{
	double V { 0. };
	for ( IntegrationPoint ip : integration_points_linear ) {
		V += ip.weight * Jacobian ( ip , vertices );
	}
	return V;
};

// PermuteNode
inline const int
Hex::PermuteNode( int node_index , int permutation_number ) {
	return permutation_to_indices(permutation_number)[node_index];
};

// PermuteCoordinates
// Takes an 8 long list of array<Point,8> and permutes according to rule
// If (bool) rule is set true, the permutation is forward (from base to use)
// else, the permutation is backward (from use to base)
inline const std::array<Euclid::Point,8>
Hex::PermuteCoordinates ( array<Euclid::Point,8> input, int permutation_number  ) {
	array<Euclid::Point,8> output;
	array<int,8> map { permutation_to_indices( permutation_number ) };
	for (size_t i = 0; i < map.size(); ++i) output[i] = input[map[i]];
	return output;
};

// double_to_point
// Converts 3 arrays of 8 array<double,8> to one array of 8 Points
// Should be templated and belong to euclid::Point
inline constexpr array<Point,8>
Hex::doubles_to_points (array<double,8> x, array<double,8> y, array<double,8> z) {
	return array<Point,8> {{
		{x[0],y[0],z[0]},{x[1],y[1],z[1]},{x[2],y[2],z[2]},{x[3],y[3],z[3]},
		{x[4],y[4],z[4]},{x[5],y[5],z[5]},{x[6],y[6],z[6]},{x[7],y[7],z[7]}
	}};
};

// edge_to_endpoints
// Map between edge index 0-11 and endpoint vertex indices 0-7
inline constexpr array<int,2>
Hex::edge_to_endpoints ( int edge_index )
{
	array<array<int,2>,12> table {{
		{0,1},{1,2},{3,2},{0,3},{4,5},{5,6},{7,6},{4,7},{0,4},{1,5},{2,6},{3,7}
	}};
	return table[edge_index];
};

// permutation_to_indices
// This list converts corner numbers 1-8
inline constexpr array<int,8>
Hex::permutation_to_indices ( int permutation ) {
	return permutation_list[permutation];
};

// local_coordinates
// Holds standard local isoparametric coordinates for a Hexahedron
inline constexpr array<Point,8>
Hex::local_coordinates {{
	{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
	{-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}
}};

// integration_points_linear
// cf. Marcin Maździarz (2010) - Unified Isoparametric 3D LagrangeFinite Elements
inline constexpr array<IntegrationPoint,8>
Hex::integration_points_linear {{
	{-0.577350269189626,-0.577350269189626, 0.577350269189626,1.},
	{ 0.577350269189626,-0.577350269189626, 0.577350269189626,1.},
	{ 0.577350269189626, 0.577350269189626, 0.577350269189626,1.},
	{-0.577350269189626, 0.577350269189626, 0.577350269189626,1.},
	{-0.577350269189626,-0.577350269189626,-0.577350269189626,1.},
	{ 0.577350269189626,-0.577350269189626,-0.577350269189626,1.},
	{ 0.577350269189626, 0.577350269189626,-0.577350269189626,1.},
	{-0.577350269189626, 0.577350269189626,-0.577350269189626,1.}
}};

// integration_points_quadratic
// cf. Marcin Maździarz (2010) - Unified Isoparametric 3D LagrangeFinite Elements
// UNUSED FOR LINEAR ELEMENT
inline constexpr array<IntegrationPoint,27>
Hex::integration_points_quadratic {{
	{-0.7745966692, -0.7745966692,-0.7745966692, 0.1714677641},
	{ 0.7745966692, -0.7745966692,-0.7745966692, 0.1714677641},
	{-0.7745966692,  0.7745966692,-0.7745966692, 0.1714677641},
	{ 0.7745966692,  0.7745966692,-0.7745966692, 0.1714677641},
	{-0.7745966692, -0.7745966692, 0.7745966692, 0.1714677641},
	{ 0.7745966692, -0.7745966692, 0.7745966692, 0.1714677641},
	{-0.7745966692,  0.7745966692, 0.7745966692, 0.1714677641},
	{ 0.7745966692,  0.7745966692, 0.7745966692, 0.1714677641},
	{ 0.          , -0.7745966692,-0.7745966692, 0.2743484225},
	{-0.7745966692,  0.          ,-0.7745966692, 0.2743484225},
	{ 0.7745966692,  0.          ,-0.7745966692, 0.2743484225},
	{ 0.          ,  0.7745966692,-0.7745966692, 0.2743484225},
	{-0.7745966692, -0.7745966692, 0.          , 0.2743484225},
	{ 0.7745966692, -0.7745966692, 0.          , 0.2743484225},
	{-0.7745966692,  0.7745966692, 0.          , 0.2743484225},
	{ 0.7745966692,  0.7745966692, 0.          , 0.2743484225},
	{ 0.          , -0.7745966692, 0.7745966692, 0.2743484225},
	{-0.7745966692,  0.          , 0.7745966692, 0.2743484225},
	{ 0.7745966692,  0.          , 0.7745966692, 0.2743484225},
	{ 0.          ,  0.7745966692, 0.7745966692, 0.2743484225},
	{ 0.          ,  0.          ,-0.7745966692, 0.438957476 },
	{ 0.          , -0.7745966692, 0.          , 0.438957476 },
	{-0.7745966692,  0.          , 0.          , 0.438957476 },
	{ 0.7745966692,  0.          , 0.          , 0.438957476 },
	{ 0.          ,  0.7745966692, 0.          , 0.438957476 },
	{ 0.          ,  0.          , 0.7745966692, 0.438957476 },
	{ 0.          ,  0.          , 0.          , 0.7023319616}
}};

// permutation_list
// This list converts contains all 24 possible solid rotations of a 3D body
// and their respective mirrors - a total of 48 possible 3D permutations.
// Generated using MATLAB
// Note for CutFEM:
// The permutation converts from THE USE CASE to THE BASE CASE.
// 36 of the 48 are used to rotate the 14 base cases onto one of the 256 case combinations.
inline constexpr array<array<int,8>,48>
Hex::permutation_list {{
	{{0,1,2,3,4,5,6,7}},  //
	{{3,2,6,7,0,1,5,4}},  // rx90
	{{4,0,3,7,5,1,2,6}},  // ry90
	{{1,2,3,0,5,6,7,4}},  // rz90
	{{7,6,5,4,3,2,1,0}},  // rx180
	{{5,4,7,6,1,0,3,2}},  // ry180
	{{2,3,0,1,6,7,4,5}},  // rz180
	{{7,3,2,6,4,0,1,5}},  // ry90 ,rx90
	{{2,6,7,3,1,5,4,0}},  // rx90 ,rz90
	{{0,3,7,4,1,2,6,5}},  // rz90 ,rx90
	{{5,1,0,4,6,2,3,7}},  // rz90 ,ry90
	{{4,5,1,0,7,6,2,3}},  // rx270
	{{1,5,6,2,0,4,7,3}},  // ry270
	{{3,0,1,2,7,4,5,6}},  // rz270
	{{6,7,3,2,5,4,0,1}},  // ry180,rx90
	{{1,0,4,5,2,3,7,6}},  // rz180,rx90
	{{3,7,4,0,2,6,5,1}},  // rx180,ry90
	{{6,2,1,5,7,3,0,4}},  // ry90 ,rx180
	{{4,7,6,5,0,3,2,1}},  // rz90 ,rx180
	{{6,5,4,7,2,1,0,3}},  // rz90 ,ry180
	{{2,1,5,6,3,0,4,7}},  // rz270,rx90
	{{7,4,0,3,6,5,1,2}},  // rx270,ry90
	{{0,4,5,1,3,7,6,2}},  // rx180,ry90 ,rx90
	{{5,6,2,1,4,7,3,0}},  // rz90 ,rx270
	{{1,0,3,2,5,4,7,6}},  // flipx
	{{3,2,1,0,7,6,5,4}},  // flipy
	{{2,3,7,6,1,0,4,5}},  // flipx,rx90
	{{5,1,2,6,4,0,3,7}},  // flipx,ry90
	{{0,3,2,1,4,7,6,5}},  // flipx,rz90
	{{0,4,7,3,1,5,6,2}},  // ry90 ,flipx
	{{2,1,0,3,6,5,4,7}},  // rz90 ,flipx
	{{4,5,6,7,0,1,2,3}},  // flipz
	{{0,1,5,4,3,2,6,7}},  // flipy,rx90
	{{7,3,0,4,6,2,1,5}},  // flipy,ry90
	{{6,7,4,5,2,3,0,1}},  // flipx,rx180
	{{6,2,3,7,5,1,0,4}},  // flipx,ry90 ,rx90
	{{4,0,1,5,7,3,2,6}},  // flipx,rz90 ,ry90
	{{1,2,6,5,0,3,7,4}},  // flipx,rz90 ,rx90
	{{7,6,2,3,4,5,1,0}},  // rx90 ,flipy
	{{3,7,6,2,0,4,5,1}},  // ry90 ,flipx,rx90
	{{3,0,4,7,2,1,5,6}},  // rz90 ,flipx,rx90
	{{4,7,3,0,5,6,2,1}},  // ry90 ,flipx,rz90
	{{1,5,4,0,2,6,7,3}},  // rz90 ,ry90 ,flipx
	{{5,6,7,4,1,2,3,0}},  // flipz,rz90
	{{5,4,0,1,6,7,3,2}},  // flipx,rx270
	{{2,6,5,1,3,7,4,0}},  // ry90 ,flipx,rx180
	{{7,4,5,6,3,0,1,2}},  // rz90 ,flipx,rx180
	{{6,5,1,2,7,4,0,3}}   // rz90 ,flipx,rx270
}};

} // namespace FEA

#endif
