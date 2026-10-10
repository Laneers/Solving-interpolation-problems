#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

double func(double x, int test_num) {
    switch (test_num) {
    case 1:
        return x * x;

    case 2:
        return 1.0 / (1.0 + x * x);

    case 3:
        return 1.0 / std::atan(1.0 + 10.0 * x * x);

    case 4: {
        double base = 4.0 * std::pow(x, 3) + 2.0 * x * x - 4.0 * x + 2.0;
        if (base < 0) {
            base = 0;
        }
        return std::pow(base, std::sqrt(2.0)) + std::asin(1.0 / (5.0 + x - x * x)) - 5.0;
    }
    default:
        std::cerr << "Error: Invalid test number (" << test_num << ").\n";
        return 0.0;
    }
}

void generateUniformGrid(double a, double b, int n, double* x) {
    double h = (b - a) / n;
    for (int i = 0; i <= n; ++i) {
        x[i] = a + i * h;
    }
}

void generateChebyshevGrid(double a, double b, int n, double* x) {
    for (int i = 0; i <= n; ++i) {
        double fact = std::cos((2.0 * i + 1.0) * M_PI / (2.0 * (n + 1)));
        x[i] = 0.5 * (a + b) + 0.5 * (b - a) * fact;
    }
}

double LagrangeInterpolation(double x_target, const double* x, const double* y, int n) {
    double result = 0.0;
    for (int i = 0; i <= n; ++i) {
        double term = y[i];
        for (int j = 0; j <= n; ++j) {
            if (j != i) {
                term *= (x_target - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

void runTestLagrange(int test_num, double a, double b, int n) {
    double* x_unif = new double[n + 1];
    double* y_unif = new double[n + 1];
    double* x_cheb = new double[n + 1];
    double* y_cheb = new double[n + 1];

    generateUniformGrid(a, b, n, x_unif);
    generateChebyshevGrid(a, b, n, x_cheb);
    
    for (int i = 0; i <= n; ++i) {
        y_unif[i] = func(x_unif[i], test_num);
        y_cheb[i] = func(x_cheb[i], test_num);
    }

    std::string filename = "test" + std::to_string(test_num) + "_n" + std::to_string(n) + "_lagr.txt";
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not create file " << filename << "\n";
    }
    else {
        outFile << std::scientific << std::setprecision(6);
        outFile << "# " << "x_val" << "\t"
            << "y_exact" << "\t"
            << "y_interp_unif" << "\t"
            << "y_interp_cheb" << "\t"
            << "err_unif" << "\t"
            << "err_cheb" << "\n";
    }

    int checkPoints = 100;
    double max_err_unif = 0.0;
    double max_err_cheb = 0.0;

    double h = (b - a) / checkPoints;
    for (int i = 0; i <= checkPoints; ++i) {
        double x_val = a + i * h;
        double y_exact = func(x_val, test_num);

        double y_interp_unif = LagrangeInterpolation(x_val, x_unif, y_unif, n);
        double y_interp_cheb = LagrangeInterpolation(x_val, x_cheb, y_cheb, n);

        double err_unif = std::abs(y_exact - y_interp_unif);
        double err_cheb = std::abs(y_exact - y_interp_cheb);

        if (err_unif > max_err_unif) max_err_unif = err_unif;
        if (err_cheb > max_err_cheb) max_err_cheb = err_cheb;

        if (outFile.is_open()) {
            outFile << x_val << "\t"
                << y_exact << "\t"
                << y_interp_unif << "\t"
                << y_interp_cheb << "\t"
                << err_unif << "\t"
                << err_cheb << "\n";
        }
    }

    if (outFile.is_open()) {
        outFile.close();
    }

    std::cout << std::setw(6) << n
        << std::setw(25) << max_err_unif
        << std::setw(25) << max_err_cheb << "\n";

    delete[] x_unif;
    delete[] y_unif;
    delete[] x_cheb;
    delete[] y_cheb;
}

int main() {
    int test_num;
	std::cout << "Input test number: "; 
    std::cin >> test_num;

    int idx = test_num - 1;
	double x_start, x_end;
    double intervals[4][2] = { {-1.0, 1.0}, {-1.0, 1.0}, {-3.0, 3.0}, {-1.0, 1.0} };

    int node_counts[] = { 5, 10, 15, 20 };

    std::cout << "Test " << test_num << " | Interval [" << intervals[idx][0] << "; " << intervals[idx][1] << "]\n";
    std::cout << std::setw(6) << "Nodes" << std::setw(25) << "Max error (uniform)" << std::setw(25) << "Max error (Chebyshev)" << "\n";

    for (int n : node_counts) {
        runTestLagrange(test_num, intervals[idx][0], intervals[idx][1], n);
    }
    std::cout << "\n";

	return 0;
}