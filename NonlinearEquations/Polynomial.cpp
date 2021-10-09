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

	Point Polynomial::findMin(Interval interval, double gaps) const
	{
		Point min{0, INT_MAX};
		for (double i = interval.first; i <= interval.second; i += gaps) {
			auto value = getValue(i);
			if (value < min.second) {
				min.first = i;
				min.second = value;
			}
		}
		return min;
	}

	Point Polynomial::findMax(Interval interval, double gaps) const
	{
		Point max{ 0, INT_MIN };
		for (double i = interval.first; i <= interval.second; i += gaps) {
			auto value = getValue(i);
			if (value > max.second) {
				max.first = i;
				max.second = value;
			}
		}
		return max;
	}

	const Polynomial::Monomials& Polynomial::getMonomials() const
	{
		return mMonomials;
	}

}