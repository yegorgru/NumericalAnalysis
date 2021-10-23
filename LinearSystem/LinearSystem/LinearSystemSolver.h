#pragma once

#include <vector>

namespace NumericalAnalysis {

	class LinearSystemSolver
	{
	public:
		static std::vector<double> GaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b);
		static std::vector<double> TridiagonalMatrix(std::vector<std::vector<double>>& A, std::vector<double>& b);
		static std::vector<std::vector<double>> Jacobi(std::vector<std::vector<double>>& A, std::vector<double>& b, double precision);
		static std::vector<std::vector<double>> UpperRelaxation(std::vector<std::vector<double>>& A, std::vector<double>& b, double precision);
	private:
		static void CheckRange(const std::vector<std::vector<double>>& A);
		static void CheckIndependence(const std::vector<std::vector<double>>& A);
		static void CheckTridiagonal(const std::vector<std::vector<double>>& A);
		static void CheckZeroOnMainDiagonal(const std::vector<std::vector<double>>& A);
	};
}

