#ifndef FEA_QUADRATURE_HARDCODED
#define FEA_QUADRATURE_HARDCODED

#include "Integration.hpp"

namespace fea {
namespace quadrature {

typedef std::vector<fea::IntegrationPoint> IntegrationRule;

namespace segment {

/* Carlos Felippa,
 * A compendium of FEM integration formulas for symbolic work,
 * Engineering Computation, Volume 21, Number 8, 2004, pages 867-890. */

const IntegrationRule felippa2 {
	{-0.57735026918962576451,1.},
	{ 0.57735026918962576451,1.}
};
const IntegrationRule felippa3 {
	{-0.77459666924148337704,0.55555555555555555556},
	{ 0.00000000000000000000,0.88888888888888888889},
	{ 0.77459666924148337704,0.55555555555555555556}
};
const IntegrationRule felippa4 {
	{-0.86113631159405257522,0.34785484513745385737},
	{-0.33998104358485626480,0.65214515486254614263},
	{ 0.33998104358485626480,0.65214515486254614263},
	{ 0.86113631159405257522,0.34785484513745385737}
};
const IntegrationRule felippa5 {
	{-0.90617984593866399280,0.23692688505618908751},
	{-0.53846931010568309104,0.47862867049936646804},
	{ 0.00000000000000000000,0.56888888888888888889},
	{ 0.53846931010568309104,0.47862867049936646804},
	{ 0.90617984593866399280,0.23692688505618908751}
};

}

} // namespace integration
} // namespace fea

#endif