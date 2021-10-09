#include "Polynomial.h"

namespace NumericalCalculus
{

	Polynomial::Polynomial(const Monomials& monomials) 
		: mMonomials(monomials)
	{

	}

	Polynomial NumericalCalculus::Polynomial::takeDerivative() const
	{
		Polynomial derivative;
		for (const auto& [degree, coef] : mMonomials) {
			if (degree != 0) {
				derivative.mMonomials[degree - 1] = degree * coef;
			}
		}
		return derivative;
	}

	const Polynomial::Monomials& Polynomial::getMonomials() const
	{
		return mMonomials;
	}

}