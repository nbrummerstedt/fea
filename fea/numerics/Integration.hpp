#ifndef FEA_INTEGRATIONPOINT
#define FEA_INTEGRATIONPOINT

#include <vector>

#include "../euclid/geometry/Point.hpp"

namespace fea {

struct IntegrationPoint {
	Euclid::Point location;
	double weight;
};
struct IntegrationRule
{
	const char name;
	const unsigned degree;
	const double tolerance;
	const std::vector<IntegrationPoint> points;
}

std::ostream& operator << (std::ostream& os, const IntegrationPoint & ip) {
	os << ip.location<<", "<< ip.weight <<"; ";
	return os;
};
std::ostream& operator << (std::ostream& os, const std::vector<IntegrationPoint> & ip) {
	for (auto p : ip) os << "\t" << p << "\n" ;
	return os;
};


}

#endif
