#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP

#include <vector>
#include <complex>

class Hamiltoniano {

private:
    double J;
    double G;
    int N;
    int dim; // Dimensión del sistema (2^N)

    // Funciones para obtener las matrices base
    std::vector<std::complex<double>> pauliX();
    std::vector<std::complex<double>> pauliZ();
    std::vector<std::complex<double>> identidad();

    // Función para calcular el producto tensorial
    std::vector<std::complex<double>> productoTensorial(
        const std::vector<std::complex<double>>& A, int filaA, int colmA,
        const std::vector<std::complex<double>>& B, int filaB, int colmB);

    // Funciones para ordenadores X y Z
    std::vector<std::complex<double>> ordenadorX(int i);
    std::vector<std::complex<double>> ordenadorZ(int i);

public:
    // Constructor personalizado
    Hamiltoniano(double j, double g, int n);

    // Construye el Hamiltoniano completo
    std::vector<std::complex<double>> construirHamiltoniano();

    //Imprime el hamiltoniano
    void imprimirMatriz(const std::vector<std::complex<double>>& matriz, int filas, int columnas) const;
};

#endif

