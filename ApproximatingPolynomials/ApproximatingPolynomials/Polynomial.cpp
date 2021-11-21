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