#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <chrono>

std::vector<double> readVector(std::ifstream& file, int size) {
    std::vector<double> vec(size);
    for (int i = 0; i < size; ++i) {
        file >> vec[i];
    }
    return vec;
}

std::vector<std::vector<double>> readMatrix(std::ifstream& file, int rows, int cols) {
    std::vector<std::vector<double>> mat(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file >> mat[i][j];
        }
    }
    return mat;
}

std::vector<double> computeResidual(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b) {
    int m = A.size(), n = A[0].size();
    std::vector<double> r(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            r[i] += A[i][j] * x[j];
        }
        r[i] = b[i] - r[i];
    }
    return r;
}

std::vector<double> solveWithCholesky(const std::vector<std::vector<double>>& T, const std::vector<double>& q) {
    std::vector<double> u(q.size(), 0.0);
    return u;
}

void solveProblem(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& c) {
    int m = A.size();
    int n = A[0].size();
    std::vector<double> x(n, 1.0);
    double gamma = 2.0 / 3.0;
    for (int iter = 0; iter < 1000; ++iter) {
        std::vector<double> r = computeResidual(A, x, b);
        std::vector<double> q(r);
        std::vector<std::vector<double>> T = A;
        std::vector<double> u = solveWithCholesky(T, q);
        std::vector<double> g(c.size(), 0.0);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                g[j] += A[i][j] * u[i];
            }
        }
        for (int j = 0; j < n; ++j) {
            g[j] = c[j] - g[j];
        }
        std::vector<double> s(g);
        for (int i = 0; i < s.size(); ++i) {
            s[i] = -s[i];
        }
        for (int i = 0; i < n; ++i) {
            x[i] += s[i];
        }
        for (int i = 0; i < m; ++i) {
            r[i] = (1 - gamma) * r[i];
        }
        double r_norm = 0.0;
        for (double val : r) {
            r_norm += val * val;
        }
        r_norm = std::sqrt(r_norm);
        if (r_norm < 1e-6) {
            break;
        }
    }
    std::cout << "x:" << std::endl;
    for (double xi : x) {
        std::cout << xi << " ";
    }
    std::cout << std::endl;
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += c[i] * x[i];
    }
    std::cout << "c*x: " << result << std::endl;
}

int main() {
    std::ifstream file("input_data.txt");
    if (!file.is_open()) {
        std::cerr << "error of open file" << std::endl;
        return -1;
    }
    std::vector<double> c = readVector(file, 9);
    std::vector<double> b = readVector(file, 6);
    std::vector<std::vector<double>> A = readMatrix(file, 6, 9);
    solveProblem(A, b, c);
    return 0;
}
