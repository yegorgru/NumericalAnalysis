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

	bool Polynomial::isDecreasing(Interval interval, double gaps) const
	{
		Polynomial d = takeDerivative();
		return d.isNegative(interval, gaps);
	}

	bool Polynomial::changeSign(Interval interval, double gaps) const
	{
		bool sign = getValue(interval.first) > 0.0;
		for (double i = interval.first + gaps; i <= interval.second; i += gaps) {
			auto value = getValue(i);
			if (value > 0 && !sign || value < 0 && sign) {
				return true;
			}
		}
		return false;
	}

	Point Polynomial::findAbsMin(Interval interval, double gaps) const
	{
		Point min{0, INT_MAX};
		for (double i = interval.first; i <= interval.second; i += gaps) {
			auto value = std::abs(getValue(i));
			if (value < min.second) {
				min.first = i;
				min.second = value;
			}
		}
		return min;
	}

	Point Polynomial::findAbsMax(Interval interval, double gaps) const
	{
		Point max{ 0, INT_MIN };
		for (double i = interval.first; i <= interval.second; i += gaps) {
			auto value = std::abs(getValue(i));
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