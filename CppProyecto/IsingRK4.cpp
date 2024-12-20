#include "hamiltoniano.hpp"
#include <vector>
#include <complex>
#include <iostream>
#include <cmath>
#include <numeric>
#include <sys/time.h>

using namespace std::complex_literals;

// Función para obtener el tiempo actual en segundos
double seconds() {
    struct timeval tmp;
    double sec;
    gettimeofday(&tmp, nullptr);
    sec = tmp.tv_sec + static_cast<double>(tmp.tv_usec) / 1000000.0;
    return sec;
}

// Función para inicializar el vector Psi
std::vector<std::complex<double>> vectorPsiInit(int n) {
    std::vector<std::complex<double>> vec(1 << n, 0.0);
    vec[0] = 1.0;
    return vec;
}

// Derivada del vector Psi
std::vector<std::complex<double>> derivadaVectorPsi(
    const std::vector<std::complex<double>>& state,
    const std::vector<std::complex<double>>& oper) {
    std::vector<std::complex<double>> result(state.size(), 0.0);
    for (size_t i = 0; i < state.size(); ++i) {
        for (size_t j = 0; j < state.size(); ++j) {
            result[i] += -1.0i * oper[i * state.size() + j] * state[j];
        }
    }
    return result;
}

// Método de Runge-Kutta 4
std::vector<std::complex<double>> rk4(
    const std::vector<std::complex<double>>& state,
    const std::vector<std::complex<double>>& oper,
    double h) {
    auto k1 = derivadaVectorPsi(state, oper);
	
	std::vector<std::complex<double>> state_k2(state.size());
	for (size_t i = 0; i < state.size(); ++i) {
		state_k2[i] = state[i] + 0.5 * h * k1[i];
    }
	auto k2 = derivadaVectorPsi(state_k2, oper);
	
	std::vector<std::complex<double>> state_k3(state.size());
	for (size_t i = 0; i < state.size(); ++i) {
		state_k3[i] = state[i] + 0.5 * h * k2[i];
    }
	auto k3 = derivadaVectorPsi(state_k3, oper);

	std::vector<std::complex<double>> state_k4(state.size());
	for (size_t i = 0; i < state.size(); ++i) {
		state_k4[i] = state[i] + h * k3[i];
    }
	auto k4 = derivadaVectorPsi(state_k4, oper);

    std::vector<std::complex<double>> result(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        result[i] = state[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    return result;
}

int main() {
    int N = 10;
    int TIMESTEPS = 1000;
    double T_MAX = 10.0;
    double h = T_MAX / TIMESTEPS;
    Hamiltoniano hamiltoniano(1.0, 1.0, N);
    auto H = hamiltoniano.construirHamiltoniano();
    int dim = 1 << N;
	
	// Solución RK4
    double start_rk4 = seconds();
    auto vecRK4 = vectorPsiInit(N);
    std::vector<double> vectorPsiRK4(TIMESTEPS, 0.0);

    for (int i = 0; i < TIMESTEPS; ++i) {
        vectorPsiRK4[i] = std::real(vecRK4[0]);
        vecRK4 = rk4(vecRK4, H, h);
    }
	
    double end_rk4 = seconds();
    std::cout << "Duración RK4: " << (end_rk4 - start_rk4) << " s" << std::endl;
    
    return 0;
}
