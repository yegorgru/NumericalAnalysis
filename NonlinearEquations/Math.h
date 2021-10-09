#pragma once

#include <map>

#include "Polynomial.h"

namespace NumericalCalculus
{
	class Math
	{
		enum class TypeOfRoot {
			MinNegative,
			MaxPositive
		};
	public:
		static double Relaxation(Polynomial polynomial, double precision, Interval interval, TypeOfRoot type);
	};
}
