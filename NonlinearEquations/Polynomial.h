#pragma once

#include <map>
#include <cmath>

namespace NumericalCalculus
{

using LowerBorder = double;
using UpperBorder = double;
using Interval = std::pair<LowerBorder, UpperBorder>;

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
public:
	const Monomials& getMonomials() const;
private:
	Monomials mMonomials;
};

} //NumericalCalculus

