import numpy as np

def derivadaVectorPsi(state, oper):
    """
    Genera la derivada del vector Psi dada en la ecuación diferencial.

    Args:
        state (np.array): Matriz que representa el estado actual.
        oper (np.array): Matriz que representa el operador (hamiltoniano).
    """
    return -1.0j * np.matmul(oper, state)

def rk4(func, oper, state, h):
    """Aplicación del método numérico Runge-Kutta de orden 4 para resolver un problema del valor inicial.

    Args:
        func (callable): Función presentada en la ecuación diferencial.
        oper (np.array): Matriz que representa al operador.
        state (np.array): Vector que representa el estado actual.
        h (float): Variable que contiene el espaciamiento temporal.

   Returns:
        (np.array): Retorna el estado después de un paso temporal.
   """
    k_1 = h*func(oper,state)
    k_2 = h*func(oper,state+((k_1)/2))
    k_3 = h*func(oper,state+((k_2)/2))
    k_4 = h*func(oper,state+k_3)
    return state + (k_1 + 2*k_2 + 2*k_3 + k_4)*(1/6)
