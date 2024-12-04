### Estudiantes:

* Santiago Arce Díaz
* Amanda Isabel Mirault Víquez

 **Coordinador**: Marlon Brenes Navarro.


 La dinámica de un estado puro para un sistema cuántico aislado se rige bajo la ecuación Schrödinger ($\hbar = 1$)
 
\begin{equation} \frac{\partial \ket{\psi(t)}}{\partial t} = -i \hat{H} \ket{\psi(t)}, \label{Ecuación 1} \end{equation}

cuya solución formal está dada por

\begin{equation} \ket{\psi(t)} = e^{-i\hat{H}(t-t_o)} \ket{\psi(t=t_o)}. \label{Ecuación } \end{equation}

Es decir, la solución involucra resolver de manera numérica la ecuación diferencial o evaluar de alguna forma la exponencial de la matriz. La idea del proyecto es evaluar la dinámica del modelo de Ising empezando de algún estado inicial.

\begin{equation} \hat{H} = -J \sum_{i=1}^N \hat{\sigma}_i ^z \hat{\sigma}_{i+1} ^z -g \sum_{i=1}^N \hat{\sigma}_i ^x, \label{Hamiltoniano} \end{equation}

donde $J$ es una escala energética que determina la interacción ferromagnética, $g$ el parámetro energético del campo transversal y $\sigma^{\alpha}_{i}$ ($\alpha = x, y, z$) son las matrices de Pauli para el spin $i$.
### Milestones:

* Con base en argumentos numéricos, evaluar cuál de las dos metodologías es mejor implementar para la solución numérica.
* Construir la matriz Hamiltoniana mediante productos tensoriales.
* Resolver el sistema utilizando algún estado inicial y visualizar su dinámica.
* Implementar la solución en Python.
* Implementar la solución en C++.
* Encontrar una forma de paralelizar el algoritmo y evaluar la aceleración
