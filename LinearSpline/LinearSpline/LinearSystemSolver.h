#pragma once

#include <vector>

namespace NumericalAnalysis {

	class LinearSystemSolver
	{
	public:
		static std::vector<double> GaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b);
	private:
		static void CheckRange(const std::vector<std::vector<double>>& A);
		static void CheckIndependence(const std::vector<std::vector<double>>& A);
	};
}

