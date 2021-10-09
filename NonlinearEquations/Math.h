#pragma once

#include <map>

#include "Polynomial.h"

namespace NumericalCalculus
{
	using LowerBorder = int;
	using UpperBorder = int;
	using Scope = std::pair<LowerBorder, UpperBorder>;
	class Math
	{
		enum class TypeOfRoot {
			MinNegative,
			MaxPositive
		};
	public:
		static double Relaxation(Polynomial polynomial, double precision, Scope scope, TypeOfRoot type);
	};
}
