#include "hamiltoniano.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

//Constructor personalizado
Hamiltoniano::Hamiltoniano(double j, double g, int n) {
    J = j;
    G = g;
    N = n;
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
std::vector<std::complex<double>> Hamiltoniano::productoTensorial(const std::vector<std::complex<double>>& A, int filaA, int colmA, const std::vector<std::complex<double>>& B, int filaB, int colmB) {
    std::vector<std::complex<double>> resultado((filaA * filaB) * (colmA * colmB), 0.0);

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
            resultado = productoTensorial(resultado, dimLocal, dimLocal, (k == i) ? pauliX() : identidad(), 2, 2);
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
            resultado = productoTensorial(resultado, dimLocal, dimLocal, (k == i || k == i + 1) ? pauliZ() : identidad(), 2, 2);
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
        for (int f = 0; f < dim; f++) {
            for (int c = 0; c < dim; c++) {
                sumZ[f * dim + c] += -J * operadorZ[f * dim + c];
            }
        }
    }

    for (int i = 1; i <= N; i++) {
        auto operadorX = ordenadorX(i);
        for (int f = 0; f < dim; f++) {
            for (int c = 0; c < dim; c++) {
                sumX[f * dim + c] += -G * operadorX[f * dim + c];
            }
        }
    }

    // Sumar ambos términos
    std::vector<std::complex<double>> hamiltoniano(dim * dim, 0);
    for (int f = 0; f < dim; f++) {
        for (int c = 0; c < dim; c++) {
            hamiltoniano[f * dim + c] = sumZ[f * dim + c] + sumX[f * dim + c];
        }
    }
    return hamiltoniano;
}

void Hamiltoniano::imprimirMatriz(const std::vector<std::complex<double>>& matriz, int filas, int columnas) const {
	for (int f = 0; f < filas; f++) {
       	  for (int c = 0; c < columnas; c++) {
            std::cout << matriz[f * columnas + c] << " "; 
        }
        std::cout << std::endl;
    }
}
