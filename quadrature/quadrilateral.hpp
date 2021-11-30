#ifndef FEA_QUADRATURE_2D_QUADRILATERAL
#define FEA_QUADRATURE_2D_QUADRILATERAL

#include "Integration.hpp"

namespace fea {
namespace integration {

struct Symmetry {
	double weight;
	double position;
	unsigned rule;
	// Symmetry rules
	// 0: origin
	// 1: axis
	// 2: diagonal
}

namespace quad {

	struct CompressedRule {
		std::string 		description;
		unsigned 		degree;
		std::vector<Symmetry> 	points;
	};
	const static CompressedRule stroud1 {
		"Hammer & Stroud", 3,{{2,0.25,sqrt(2./3.)}}
	};
	const static CompressedRule stroud2 {
		"Hammer & Stroud", 5, {
			{0,1.,16./81.},
			{1,sqrt(3./5.),10/81.},
			{2,sqrt(3./5.),25/(4*81.)}}
	};
}

} // namespace integration
} // namespace fea

#endif
