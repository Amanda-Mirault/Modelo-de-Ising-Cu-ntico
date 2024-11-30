#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP

#include <vector>
#include <complex>

class Hamiltoniano {

private:
    double J;
    double g;
    int N;
    int dim; // Dimensión del sistema (2^N)

    // Funciones privadas para obtener matrices base
    std::vector<std::complex<double>> pauliX();
    std::vector<std::complex<double>> pauliZ();
    std::vector<std::complex<double>> identidad();

    // Función para calcular el producto tensorial con índices lineales
    std::vector<std::complex<double>> productoTensorial(
        const std::vector<std::complex<double>>& A, int filaA, int colmA,
        const std::vector<std::complex<double>>& B, int filaB, int colmB
    );

    // Funciones para operadores X y Z con índices lineales
    std::vector<std::complex<double>> ordenadorX(int i);
    std::vector<std::complex<double>> ordenadorZ(int i);

public:
    // Constructor personalizado
    Hamiltoniano(double J, double g, int N);

    // Construye el Hamiltoniano completo
    std::vector<std::complex<double>> construirHamiltoniano();
};

#endif

