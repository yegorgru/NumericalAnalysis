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
		static Polynomial multiply(const Polynomial& lhs, const Polynomial& rhs);
		static Polynomial add(const Polynomial& lhs, const Polynomial& rhs);
	public:
		Polynomial() = default;
		Polynomial(const Monomials& monomials);
	public:
		double getValue(double x) const;
		const Monomials& getMonomials() const;
		bool equal(const Polynomial& rhs, double precision) const;
	private:
		Monomials mMonomials;
	};

} //NumericalAnalysis