#pragma once

#include <map>
#include <cmath>

namespace NumericalAnalysis
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
	private:
		Monomials mMonomials;
	};

} //NumericalAnalysis