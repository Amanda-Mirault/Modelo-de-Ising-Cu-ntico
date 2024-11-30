#include "hamiltoniano.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

//Constructor personalizado
Hamiltoniano::Hamiltoniano(double J, double g, int N) {
    J = J;
    g = g;
    N = N;
    dim = 1 << N;
}

// Matrices de Pauli y matriz identidad
std::vector<std::complex<double>> Hamiltoniano::pauliX() {
    return {0, 1, 1, 0};
}

std::vector<std::complex<double>> Hamiltoniano::pauliZ() {
    return {1, 0, 0, -1};
}

std::vector<std::complex<double>> Hamiltoniano::identidad() {
    return {1, 0, 0, 1};
}

//Producto tensorial
std::vector<double> productoTensorial(const std::vector<double>& A, int filaA, int colmA, const std::vector<double>& B, int filaB, int colmB) {
    std::vector<double> resultado((filaA * filaB) * (colmA * colmB), 0.0);

    for (int i = 0; i < filaA; ++i) {
        for (int j = 0; j < colmA; ++j) {
            for (int k = 0; k < filaB; ++k) {
                for (int l = 0; l < colmB; ++l) {
                    resultado[(i * filaB + k) * (colmA * colmB) + (j * colmB + l)] = A[i * colmA + j] * B[k * colmB + l];
                }
            }
        }
    }

    return resultado;
}

//Ordenador para X
std::vector<std::complex<double>> Hamiltoniano::ordenadorX(int i) {
	
    std::vector<std::complex<double>> resultado = identidad();
    int dimLocal = 2;

    for (int k = 1; k <= N; k++) {
        if (k == 1) {
            resultado = (k == i) ? pauliX() : identidad();
        } else {
            resultado = productoTensorial(resultado, (k == i) ? pauliX() : identidad(), dimLocal, 2);
            dimLocal *= 2;
        }
    }
    return resultado;
}

//Ordenador para Z
std::vector<std::complex<double>> Hamiltoniano::ordenadorZ(int i) {
    std::vector<std::complex<double>> resultado = identidad();
    int dimLocal = 2;

    for (int k = 1; k <= N; k++) {
        if (k == 1) {
            resultado = (k == i) ? pauliZ() : identidad();
        } else {
            resultado = productoTensorial(resultado, (k == i || k == i + 1) ? pauliZ() : identidad(), dimLocal, 2);
            dimLocal *= 2;
        }
    }
    return resultado;
}

//Se construye el Hamiltoniano
std::vector<std::complex<double>> Hamiltoniano::construirHamiltoniano() {
	std::vector<std::complex<double>> sumZ(dim * dim, 0);
	std::vector<std::complex<double>> sumX(dim * dim, 0);

    for (int i = 1; i < N; i++) {
        auto operadorZ = ordenadorZ(i);
        for (int r = 0; r < dim; r++) {
            for (int c = 0; c < dim; c++) {
                sumZ[r * dim + c] += -J * operadorZ[r * dim + c];
            }
        }
    }

    for (int i = 1; i <= N; i++) {
        auto operadorX = ordenadorX(i);
        for (int r = 0; r < dim; r++) {
            for (int c = 0; c < dim; c++) {
                sumX[r * dim + c] += -g * operadorX[r * dim + c];
            }
        }
    }

    // Sumar ambos términos
    vector<complex<double>> hamiltoniano(dim * dim, 0);
    for (int r = 0; r < dim; r++) {
        for (int c = 0; c < dim; c++) {
            hamiltoniano[r * dim + c] = sumZ[r * dim + c] + sumX[r * dim + c];
        }
    }
    return hamiltoniano;
}
