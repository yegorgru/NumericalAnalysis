#pragma once

#include <map>
#include <cmath>

namespace NumericalCalculus
{

using Interval = std::pair<double, double>;
using Point = std::pair<double, double>;

class Polynomial
{
public:
	using Coef = double;
	using Degree = int;
	using Monomials = std::map<Degree, Coef>;
public:
	Polynomial() = default;
	Polynomial(const Monomials& monomials);
public:
	double getValue(double x) const;
	Polynomial takeDerivative() const;
	bool isPositive(Interval interval, double gaps) const;
	bool isNegative(Interval interval, double gaps) const;
	bool isIncreasing(Interval interval, double gaps) const;
	bool isDecreasing(Interval interval, double gaps) const;
	bool changeSign(Interval interval, double gaps) const;
	Point findAbsMin(Interval interval, double gaps) const;
	Point findAbsMax(Interval interval, double gaps) const;
public:
	const Monomials& getMonomials() const;
private:
	Monomials mMonomials;
};

} //NumericalCalculus

