import numpy as np

class Hamiltoniano:
  def __init__(self, J, g, N):
    self.J = J
    self.g = g
    self.N = N

  def pauliX(self):
    return np.array([[0,1],[1,0]])

  def pauliZ(self):
      return np.array([[1,0],[0,-1]])

  def identidad(self):
      return np.array([[1,0],[0,1]])

  def productoTensorial(self, matrizA, matrizB):
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

  def ordenadorX(self, i):
    result = self.identidad()

    for k in range(1, self.N + 1):
        if k == 1:
            if k == i:
                result = self.pauliX()
            else:
                result = self.identidad()
        else:
            if k == i:
                result = self.productoTensorial(result, self.pauliX())
            else:
                result = self.productoTensorial(result, self.identidad())
    return result

  def ordenadorZ(self, i):
    result = self.identidad()

    for k in range(1, self.N + 1):
        if k == 1:
            if k == i:
                result = self.pauliZ()
            else:
                result = self.identidad()
        elif k == i:
            result = self.productoTensorial(result, self.pauliZ())
        elif k == i+1:
            result = self.productoTensorial(result, self.pauliZ())
        else:
            result = self.productoTensorial(result, self.identidad())
    return result

  def construirHamiltoniano(self):
    sumZ = np.zeros((2**self.N, 2**self.N))
    for i in range(1, self.N):
        sumZ += self.ordenadorZ(i)
    sumZ *= -self.J

    sumX = np.zeros((2**self.N, 2**self.N))
    for i in range(1, self.N + 1):
        sumX += self.ordenadorX(i)
    sumX *= -self.g

    return sumZ + sumX

H = Hamiltoniano(1, 1, 2)
print(H.construirHamiltoniano())
