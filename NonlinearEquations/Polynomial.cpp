#include "Polynomial.h"

namespace NumericalCalculus
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

	bool Polynomial::isPositive(Interval interval, double gaps) const
	{
		for (double i = interval.first; i <= interval.second; i += gaps) {
			if (getValue(i) < 0) {
				return false;
			}
		}
		return true;
	}

	bool Polynomial::isNegative(Interval interval, double gaps) const
	{
		for (double i = interval.first; i <= interval.second; i += gaps) {
			if (getValue(i) > 0) {
				return false;
			}
		}
		return true;
	}

	bool Polynomial::isIncreasing(Interval interval, double gaps) const
	{
		Polynomial d = takeDerivative();
		return d.isPositive(interval, gaps);
	}

	const Polynomial::Monomials& Polynomial::getMonomials() const
	{
		return mMonomials;
	}

}