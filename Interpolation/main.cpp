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
        double fact = - std::cos((2.0 * i + 1.0) * M_PI / (2.0 * (n + 1)));     //in ascending order
        x[i] = 0.5 * (a + b) + 0.5 * (b - a) * fact;
    }
}

double LagrangeInterpolation(double x_target, const double* x, const double* y, int n) {
    //Basic form
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

double SplineInterpolation(double x_target, const double* x, int n, const double* a, const double* b, const double* c, const double* d) {
    int idx = 1;
    for (int i = 1; i <= n; ++i) {
        if (x_target <= x[i]) {
            idx = i;
            break;
        }
        if (i == n) idx = n;
    }
    double x_i = x_target - x[idx - 1];

    double s = a[idx] + b[idx] * x_i + c[idx] * x_i * x_i + d[idx] * x_i * x_i * x_i;

    return s;
}

void Progon(const double* h, const double* y, int n, double* c) {
    c[1] = 0;
    c[n + 1] = 0;

    double* alpha = new double[n + 1];
    double* beta = new double[n + 1];

    alpha[1] = 0.0;
    beta[1] = 0.0;

    for (int i = 2; i <= n; ++i) {
        double A = h[i - 1];
        double B = 2.0 * (h[i - 1] + h[i]);
        double C = h[i];
        double F = 3.0 * ((y[i] - y[i - 1]) / h[i] - (y[i - 1] - y[i - 2]) / h[i - 1]);

        double denom = A * alpha[i - 1] + B;
        alpha[i] = -C / denom;
        beta[i] = (F - A * beta[i - 1]) / denom;
    }

    for (int i = n; i >= 2; --i) {
        c[i] = alpha[i] * c[i + 1] + beta[i];
    }

    delete[] alpha;
    delete[] beta;
}

void SplineCoeffs(const double* x, const double* y, int n, double* a, double* b, double* c, double* d) {
    double* h = new double[n + 1];
    for (int i = 1; i < n + 1; ++i) {
        h[i] = x[i] - x[i - 1];
        a[i] = y[i - 1];
    }
    Progon(h, y, n, c);
    for (int i = 1; i < n + 1; ++i) {
        b[i] = (y[i] - y[i - 1]) / h[i] - c[i] * h[i] - (c[i + 1] - c[i]) * h[i]/ 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    delete h;
}

void runTest(int test_num, double a, double b, int n) {
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

    //Spline coeffs
    double* a_unif = new double[n + 1];
    double* b_unif = new double[n + 1];
    double* c_unif = new double[n + 2];
    double* d_unif = new double[n + 1];
    SplineCoeffs(x_unif, y_unif, n, a_unif, b_unif, c_unif, d_unif);

    double* a_cheb = new double[n + 1];
    double* b_cheb = new double[n + 1];
    double* c_cheb = new double[n + 2];
    double* d_cheb = new double[n + 1];
    SplineCoeffs(x_cheb, y_cheb, n, a_cheb, b_cheb, c_cheb, d_cheb);

    std::string filename_lagr = "test" + std::to_string(test_num) + "_n" + std::to_string(n) + "_lagr.txt";
    std::ofstream out_lagr(filename_lagr);

    std::string filename_spline = "test" + std::to_string(test_num) + "_n" + std::to_string(n) + "_spline.txt";
    std::ofstream out_spline(filename_spline);

    out_lagr << std::scientific << std::setprecision(6);
    out_lagr << std::setw(14) << "x_val" << "\t"
        << std::setw(14) << "y_exact" << "\t"
        << std::setw(14) << "y_interp_unif" << "\t"
        << std::setw(14) << "y_interp_cheb" << "\t"
        << std::setw(14) << "err_unif" << "\t"
        << std::setw(14) << "err_cheb" << "\n";

    out_spline << std::scientific << std::setprecision(6);
    out_spline << std::setw(14) << "x_val" << "\t"
        << std::setw(14) << "y_exact" << "\t"
        << std::setw(14) << "y_spline_unif" << "\t"
        << std::setw(14) << "y_spline_cheb" << "\t"
        << std::setw(14) << "err_unif" << "\t"
        << std::setw(14) << "err_cheb" << "\n";

    int n_checks = 100;
    double max_err_lagr_unif = 0.0;
    double max_err_lagr_cheb = 0.0;
    double max_err_spline_unif = 0.0;
    double max_err_spline_cheb = 0.0;

    double h = (b - a) / n_checks;
    for (int i = 0; i <= n_checks; ++i) {
        double x_val = a + i * h;
        double y_exact = func(x_val, test_num);

        double y_lagr_unif = LagrangeInterpolation(x_val, x_unif, y_unif, n);
        double y_lagr_cheb = LagrangeInterpolation(x_val, x_cheb, y_cheb, n);
        double err_lagr_unif = std::abs(y_exact - y_lagr_unif);
        double err_lagr_cheb = std::abs(y_exact - y_lagr_cheb);

        if (err_lagr_unif > max_err_lagr_unif) max_err_lagr_unif = err_lagr_unif;
        if (err_lagr_cheb > max_err_lagr_cheb) max_err_lagr_cheb = err_lagr_cheb;

        out_lagr << x_val << "\t"
            << y_exact << "\t"
            << y_lagr_unif << "\t"
            << y_lagr_cheb << "\t"
            << err_lagr_unif << "\t"
            << err_lagr_cheb << "\n";

        double y_spline_unif = SplineInterpolation(x_val, x_unif, n, a_unif, b_unif, c_unif, d_unif);
        double y_spline_cheb = SplineInterpolation(x_val, x_cheb, n, a_cheb, b_cheb, c_cheb, d_cheb);
        double err_spline_unif = std::abs(y_exact - y_spline_unif);
        double err_spline_cheb = std::abs(y_exact - y_spline_cheb);

        if (err_spline_unif > max_err_spline_unif) max_err_spline_unif = err_spline_unif;
        if (err_spline_cheb > max_err_spline_cheb) max_err_spline_cheb = err_spline_cheb;

        out_spline << x_val << "\t"
            << y_exact << "\t"
            << y_spline_unif << "\t"
            << y_spline_cheb << "\t"
            << err_spline_unif << "\t"
            << err_spline_cheb << "\n";
    }

    out_lagr.close();
    out_spline.close();

    std::cout << std::setw(6) << n
        << std::setw(20) << max_err_lagr_unif
        << std::setw(20) << max_err_lagr_cheb
        << std::setw(20) << max_err_spline_unif
        << std::setw(20) << max_err_spline_cheb << "\n";

    delete[] x_unif;
    delete[] y_unif;
    delete[] x_cheb;
    delete[] y_cheb;
}

int findNNodes(int test_num, double a, double b, double eps, int max_n = 200) {
    int n_checks = 1000;
    double h = (b - a) / n_checks;

    for (int n = 1; n <= max_n; ++n) {
        double* x_cheb = new double[n + 1];
        double* y_cheb = new double[n + 1];

        generateChebyshevGrid(a, b, n, x_cheb);
        for (int i = 0; i <= n; ++i) {
            y_cheb[i] = func(x_cheb[i], test_num);
        }

        double max_err = 0.0;
        for (int i = 0; i <= n_checks; ++i) {
            double x_val = a + i * h;
            double y_exact = func(x_val, test_num);
            double y_interp = LagrangeInterpolation(x_val, x_cheb, y_cheb, n);

            double err = std::abs(y_exact - y_interp);
            if (err > max_err) {
                max_err = err;
            }
        }

        delete[] x_cheb;
        delete[] y_cheb;

        if (max_err <= eps) {
            return n + 1;
        }
    }

    return -1;
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
    std::cout << std::setw(6) << "Nodes" << std::setw(20) << "Lagr (unif)" << std::setw(20) << "Lagr (cheb)" << std::setw(20) << "Spline (unif)" << std::setw(20) << "Spline (cheb)" << "\n";

    for (int n : node_counts) {
        runTest(test_num, intervals[idx][0], intervals[idx][1], n);
    }
    std::cout << "\n";

    double eps_list[] = { 0.01, 0.0001 };
    std::cout << "Required Chebyshev nodes for Lagrange interpolation:\n";
    for (double eps : eps_list) {
        int nodes = findNNodes(test_num, intervals[idx][0], intervals[idx][1], eps);
        std::cout << "  eps = " << std::defaultfloat << eps << " -> ";
        if (nodes != -1) {
            std::cout << nodes << " nodes (n = " << nodes - 1 << ")\n";
        }
        else {
            std::cout << "more than 200 nodes required\n";
        }
    }

	return 0;
}