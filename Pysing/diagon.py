import munpy as np

def diagonalizacion(t, t_0, vec, eigVals, eigVecs):
    """Solución directa de la ecuación diferencial por medio de digonalización.

    Args:
        t (float): Variable que representa el tiempo en el estado actual.
        t_0 (float): Variable que representa el tiempo en el estado inicial.
        vec (np.array): Vector en el estado inicial.
        eigVals (np.array): Vector que contiene los autovalores.
        eigVecs (np.array): Matriz que contiene los autovectores respectivos.

    Returns:
        (np.array): Retorna el estado del vector en el tiempo t.
    """
    eigVals = -1.0j * eigVals * (t - t_0)
    expVals = np.exp(eigVals)
    expD = np.diag(expVals)

    psi_d = np.dot(eigVecs.T, vec)
    psi_t_d = np.dot(expD, psi_d)

    psi_t = np.dot(eigVecs, psi_t_d)

    return psi_t
