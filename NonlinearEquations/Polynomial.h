#pragma once

#include <map>

namespace NumericalCalculus
{

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
	Polynomial takeDerivative() const;
public:
	const Monomials& getMonomials() const;
private:
	Monomials mMonomials;
};

} //NumericalCalculus

