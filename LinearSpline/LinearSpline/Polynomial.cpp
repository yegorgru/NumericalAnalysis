#include "Polynomial.h"

namespace NumericalAnalysis
{

	Polynomial::Polynomial(const Monomials& monomials)
		: mMonomials(monomials)
	{

	}

	double Polynomial::getValue(double x) const
	{
		double value = 0.0;
		for (const auto& [degree, coef] : mMonomials) {
			value += std::pow(x, degree) * coef;
		}
		return value;
	}
}