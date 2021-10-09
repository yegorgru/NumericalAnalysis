#pragma once

#include <map>
#include <vector>
#include <iostream>

#include "Polynomial.h"

namespace NumericalCalculus
{
	class Math
	{
	public:
		static std::vector<Point> Relaxation(Polynomial polynomial, double precision, Interval interval);
	};
}
