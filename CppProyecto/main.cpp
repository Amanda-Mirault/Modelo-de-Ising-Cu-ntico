#include "hamiltoniano.hpp"
#include <iostream>
#include <vector>
#include <complex>

int main() {
	int N = 2;
	Hamiltoniano Ham(1.0, 1.0, N);
	std::vector<std::complex<double>> H = Ham.construirHamiltoniano();
	std::cout << "Hamiltoniano:" << std::endl;

	int dim = 1 << N;
	Ham.imprimirMatriz(H, dim, dim);

	return 0;
}
