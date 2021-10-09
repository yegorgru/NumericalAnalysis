#pragma once

#include <map>

#include "Polynomial.h"

namespace NumericalCalculus
{
	class Math
	{
	public:
		static Point Relaxation(Polynomial polynomial, double precision, Interval interval);
	};
}
