import numpy as np

class Hamiltoniano:
    """
    Generador de una matriz Hamiltoniana a partir de las matrices de Pauli X y Z.

    Attributes:
        J (float): Escala energética que determina la interacción ferromagnética. 
        g (float): Parámetro energéticos del campo transversal.
        N (int): Determina la cantidad de spines.

    Methods:
        pauliX(): Matriz X de Pauli.
        pauliZ(): Matriz Z de Pauli.
        identidad(): Matriz identidad.
        productoTensorial(matrizA, matrizB): Calcula el producto tensorial entre dos matrices.
        ordenadorX(i): Implementa el orden correcto de los productos tensoriales de las matrices de Pauli X y la identidad.
        ordenadorZ(i): Implementa el orden correcto de los productos tensoriales de las matrices de Pauli Z y la identidad.
        construirHamiltoniano(): Construye el hamiltoniano final.
    """
    def __(self, J, g, N):
        """
        Inicializa un objeto de tipo Hamiltoniano.

        Args:
            J (float): Escala energética que determina la interacción ferromagnética. 
            g (float): Parámetro energéticos del campo transversal.
            N (int): Determina la cantidad de spines.
        """
        self.J = J
        self.g = g
        self.N = N

    def pauliX(self):
        """
        Genera la matriz X de Pauli

        Returns:
            (np.array): Retorna la matriz X de Pauli.
        """
        return np.array([[0,1],[1,0]])

    def pauliZ(self):     
        """
        Genera la matriz Z de Pauli

        Returns:
            (np.array): Retorna la matriz Z de Pauli.
        """
        return np.array([[1,0],[0,-1]])

    def identidad(self):
        """
        Genera la matriz identidad.

        Returns:
            (np.array): Retorna la matriz identidad.
        """
        return np.array([[1,0],[0,1]])

    def productoTensorial(self, matrizA, matrizB):
        """
        Realiza el producto tensorial entre la matriz A y la matriz B.
        Example:
            >>> import numpy as np
            >>> productoTensorial(np.array([[1, 2], [3, 4]]), np.array([[0, 5], [6, 7]]))
            [[ 0.  5.  0. 10.]
             [ 6.  7. 12. 14.]
             [ 0. 15.  0. 20.]
             [18. 21. 24. 28.]]

        Args:
            matrizA (np.array): Matriz A
            matrizB (np.array): Matriz B

        Returns:
            (np.array): Retorna el producto tensorial entre la matriz A y B.
        """
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
        """
        Implementa el orden correcto de los productos tensoriales de las matrices de Pauli X y la identidad.
        
        Args:
            i (int): Indica la posición en el producto tensorial de la matriz X de Pauli.

        Returns:
            (np.array): Retorna el producto tensorial entre la matriz X de Pauli y la identidad según la posición i.
        """
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
        """
        Implementa el orden correcto de los productos tensoriales de las matrices de Pauli Z y la identidad.

        Args:
            i (int): Indica la posición en el producto tensorial de las dos matrices Z de Pauli.

        Returns:
            (np.array): Retorna el producto tensorial entre las matrices Z de Pauli y la identidad según la posición i.
        """

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
        """
        Construye el hamiltoniano para el modelo de Ising.

        Example:
            >>> Ham = Hamiltoniano(1, 1, 2)
            >>> Ham.construirHamiltoniano()
            [[-1. -1. -1. -0.]
             [-1.  1. -0. -1.]
             [-1. -0.  1. -1.]
             [-0. -1. -1. -1.]]

        Returns:
            (np.array): Retorna el hamiltoniano para el modelo de Ising de dimensión 2^N * 2^N.
        """
        sumZ = np.zeros((2**self.N, 2**self.N))
        for i in range(1, self.N):
            sumZ += self.ordenadorZ(i)
        sumZ *= -self.J

        sumX = np.zeros((2**self.N, 2**self.N))
        for i in range(1, self.N + 1):
            sumX += self.ordenadorX(i)
        sumX *= -self.g

        return sumZ + sumX
