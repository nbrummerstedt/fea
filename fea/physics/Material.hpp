#ifndef FEA_MATERIAL
#define FEA_MATERIAL

#include <Eigen/Dense>
#include "../euclid/Geometry"

namespace FEA {
class Material {
	double _young;
	double _poisson;
	public:
	Material(double E=1.,double nu=0.3) : _young(E), _poisson(nu) {};
	const Eigen::Matrix<double,6,6> ConstitutiveLinearIsotropicElastic ();
};

inline const Eigen::Matrix<double,6,6>
Material::ConstitutiveLinearIsotropicElastic () {
	double l = ( _young * _poisson ) / ((1. + _poisson ) * (1. - 2. * _poisson ));
	double g =   _young / (2. * (1. + _poisson ));
	double G = 2. * g;
	Eigen::Matrix<double,6,6> C;
	C << 	l+G, 	l, 		l, 		0., 	0., 	0.,
			l, 		l+G, 	l, 		0., 	0., 	0.,
			l, 		l, 		l+G, 	0., 	0., 	0.,
			0., 	0., 	0., 	g, 		0., 	0.,
			0., 	0., 	0., 	0., 	g, 		0.,
			0., 	0., 	0., 	0., 	0., 	g;
	return C;
};

}

#endif
