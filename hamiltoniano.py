import numpy as np

#matricitas necesarias
def pauliX():
    return np.array([[0,1],[1,0]])

def pauliZ():
    return np.array([[1,0],[0,-1]])

def identidad():
    return np.array([[1,0],[0,1]])

#funcion para calcular producto tensorial
def productoTensorial(matrizA, matrizB):
    m = len(matrizA)
    n = len(matrizA[0])
    p = len(matrizB)
    q = len(matrizB[0])

    result = np.zeros((m * p, n * q))

    for i in range(m):
        for j in range(n):
            for k in range(p):
                for l in range(q):
                    result[i * p + k][j * q + l] = matrizA[i][j] * matrizB[k][l]

    return result

#ordenar productos tensoriales para x
def ordenadorX(N, i, sigma, base):
    result = base()

    for k in range(1, N+1):
        if k == 1:
            if k == i:
                result = sigma()
            else:
                result = base()
        else:
            if k == i:
                result = productoTensorial(result, sigma())
            else:
                result = productoTensorial(result, base())
    return result

#ordenar productos tensoriales para z
def ordenadorZ(N, i, sigma, base):
    result = base()

    for k in range(1, N+1):
        if k == 1:
            if k == i:
                result = sigma()
            else:
                result = base()
        elif k == i:
            result = productoTensorial(result, sigma())
        elif k == i+1:
            result = productoTensorial(result, sigma())
        else:
            result = productoTensorial(result, base())
    return result

#construir hamiltoniano
def hamiltoniano(N, J, g, sigmaZ, sigmaX, base):
    sumZ = np.zeros((2**N, 2**N))
    for i in range(1,N):
        sumZ += ordenadorZ(N, i, sigmaZ, base)
    sumZ *= -J

    sumX = np.zeros((2**N, 2**N))
    for i in range(1, N+1):
        sumX += ordenadorX(N, i, sigmaX, base)
    sumX *= -g

    return sumZ + sumX

H = hamiltoniano(2, 1, 1, pauliZ, pauliX, identidad)

print("Hamiltoniano: ", H)
