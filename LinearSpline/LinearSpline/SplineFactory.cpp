#include "SplineFactory.h"

#include "LinearSystemSolver.h"

#include <stdexcept>
#include <iostream>

namespace NumericalAnalysis {
	
	SplineFactory::Spline SplineFactory::createLinearSpline(const Polynomial& function, range_t range, uint32_t nodes, bool debug)
	{
		if (nodes < 2) {
			throw std::runtime_error("Incorrect number of nodes");
		}
		double step = (range.second - range.first) / (nodes - 1);
		double x0 = range.first;
		double x1 = x0 + step;
		Spline result;
		while (x1 <= range.second) {
			std::vector<std::vector<double>> A{
				{x0, 1},
				{x1, 1}
			};
			double f0 = function.getValue(x0);
			double f1 = function.getValue(x1);
			std::vector<double> b{ f0, f1 };
			if (debug) {
				std::cout << "\nRange: (" << x0 << ", " << x1 << "):\n"
					<< "System:\n"
					<< x0 << "a + b = " << f0 << "\n"
					<< x1 << "a + b = " << f1 << std::endl;
			}
			auto coefficients = LinearSystemSolver::GaussianElimination(A, b);
			if (debug) {
				std::cout << "Coefficients: a = " << coefficients.at(0) << ", b = " << coefficients.at(1) << "\n"
					<< "Function: y = " << coefficients.at(0) << "x + " << coefficients.at(1) << std::endl;
			}
			result.emplace(std::make_pair(x0, x1), Polynomial ({
				{1, coefficients.at(0)},
				{0, coefficients.at(1)}
			}));
			x0 += step;
			x1 += step;
		}
		return result;
	}

} //namespace NumericalAnalysis