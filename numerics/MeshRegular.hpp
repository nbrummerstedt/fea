#ifndef FEA_MESHREGULAR_H
#define FEA_MESHREGULAR_H

#include <cmath>
#include <iostream>

namespace fea {

class MeshRegular {

	double xmin, xmax, ymin, ymax, zmin, zmax;
	unsigned resolution;

	public:
	MeshRegular ( double x1, double x2, double y1, double y2, double z1, double z2, unsigned R )
		: xmin(x1),xmax(x2),ymin(y1),ymax(y2),zmin(z1),zmax(z2),resolution(R) {};

	unsigned nEx() { return ( xmax-xmin >= ymax-ymin && xmax-xmin >= zmax-zmin) ? R : ( ymax-ymin >= xmax-xmin && ymax-ymin >= zmax-zmin) ? ceil( (double) R * (xmax-xmin)/(ymax-ymin) ) : ( zmax-zmin >= xmax-xmin && zmax-zmin >= ymax-ymin) ? ceil( (double) R * (xmax-xmin)/(zmax-zmin) ) : 0 ); }
	unsigned nEy() { return ( xmax-xmin >= ymax-ymin && xmax-xmin >= zmax-zmin) ? ceil( (double) R * (ymax-ymin)/(xmax-xmin) ) : ( ymax-ymin >= xmax-xmin && ymax-ymin >= zmax-zmin) ? R : ( zmax-zmin >= xmax-xmin && zmax-zmin >= ymax-ymin) ? ceil( (double) R * (ymax-ymin)/(zmax-zmin) ) : 0 ); }
	unsigned nEz() { return	( xmax-xmin >= ymax-ymin && xmax-xmin >= zmax-zmin) ? ceil( (double) R * (zmax-zmin)/(xmax-xmin) ) : ( ymax-ymin >= xmax-xmin && ymax-ymin >= zmax-zmin) ? ceil( (double) R * (zmax-zmin)/(ymax-ymin) ) : ( zmax-zmin >= xmax-xmin && zmax-zmin >= ymax-ymin) ? R : 0 ); }

	unsigned nVx() { return nEx()+1; };
	unsigned nVy() { return nEy()+1; };
	unsigned nVz() { return nEz()+1; };
	unsigned VertexCount()  { return nVx() * nVy() * nVz(); }
	unsigned ElementCount() { return nEx*nEy*nEz; }

	double dx() { return xmax-xmin; };
	double dy() { return ymax-ymin; };
	double dz() { return zmax-zmin; };

	Euclid::Point
	Vertex( unsigned vertex_index ) {
		// Quit early if input is wrong
		if ( vertex_index + 1 > VertexCount() ) { throw("Index out of bounds."); };
		// Find z-slice, y-line, and x-number
		unsigned nVZ { (nVx() * nVy()) };
		unsigned vertex_count_line { nVx() };
		unsigned zn { ( vertex_index / nVZ ) };
		unsigned yn { ( vertex_index % nVZ ) / vertex_count_line };
		unsigned xn { ( vertex_index % nVZ ) % vertex_count_line };
		// Convert to coordinates
		double z { zmin + ( (double) zn / (double) (nEz) ) * (zmax-zmin) };
		double y { ymin + ( (double) yn / (double) (nEy) ) * (zmax-zmin) };
		double x { xmin + ( (double) xn / (double) (nEx) ) * (zmax-zmin) };
		return Euclid::Point {x,y,z};
	};

};

}

#endif
