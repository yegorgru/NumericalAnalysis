#pragma once

#include "Polynomial.h"

#include <map>
#include <utility>

namespace NumericalAnalysis {
	using range_t = std::pair<double, double>;

	class SplineFactory
	{
	public:
		using Spline = std::map<range_t, Polynomial>;
		static Spline createLinearSpline(const Polynomial& function, range_t range, uint32_t nodes, bool debug);
	};
} //namespace NumericalAnalysis
