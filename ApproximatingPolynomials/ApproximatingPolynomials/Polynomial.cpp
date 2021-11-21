#include "Polynomial.h"

namespace NumericalAnalysis
{

	Polynomial Polynomial::multiply(const Polynomial& lhs, const Polynomial& rhs)
	{
		Polynomial::Monomials result;
		for (const auto& lMonomial : lhs.getMonomials()) {
			for (const auto& rMonomial : rhs.getMonomials()) {
				result[lMonomial.first + rMonomial.first] += lMonomial.second * rMonomial.second;
			}
		}
		for (auto it = result.begin(); it != result.end(); ++it) {
			if (it->second == 0) {
				it = result.erase(it);
			}
		}
		return Polynomial(result);
	}

	Polynomial Polynomial::add(const Polynomial& lhs, const Polynomial& rhs)
	{
		Polynomial::Monomials result = lhs.getMonomials();
		for (const auto& rMonomial : rhs.getMonomials()) {
			result[rMonomial.first] += rMonomial.second;
		}
		for (auto it = result.begin(); it != result.end(); ++it) {
			if (it->second == 0) {
				it = result.erase(it);
			}
		}
		return Polynomial(result);
	}

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

	Polynomial Polynomial::takeDerivative() const
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
		Point min{ 0, INT_MAX };
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

	bool Polynomial::equal(const Polynomial& rhs, double precision) const
	{
		const auto& rMonomials = rhs.getMonomials();
		auto it1 = mMonomials.begin();
		auto it2 = rMonomials.begin();
		for (; it1 != mMonomials.end() && it2 != rMonomials.end(); ++it1, ++it2) {
			if (std::abs(it1->second - it2->second) > precision || it1->first != it2->first) {
				return false;
			}
		}
		if (it1 != mMonomials.end() || it2 != rMonomials.end()) {
			return false;
		}
		return true;
	}
}