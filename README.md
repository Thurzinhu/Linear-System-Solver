# Métodos de Aproximação de Sistemas Lineares

Esta implementação inclui cinco classes em Python que oferecem diferentes métodos para a aproximação de soluções de sistemas lineares. Estas classes foram desenvolvidas por Arthur Andrade D'Olival em outubro de 2023.

## Classes

1. `Método`:
   - Atributos:
     - `coeficientes`: Matriz de coeficientes para o sistema linear.
     - `resultado`: Vetor de resultados para o sistema linear.

2. `Jacobi`:
   - Herda de `Método`.
   - Implementa o método de Jacobi para a aproximação de soluções de sistemas lineares.

3. `Seidel`:
   - Herda de `Jacobi`.
   - Implementa o método de Gauss-Seidel, uma variação do método de Jacobi.

4. `Gauss`:
   - Herda de `Método`.
   - Implementa o método de eliminação gaussiana seguido de substituição retroativa para encontrar soluções de sistemas lineares.

5. `LU`:
   - Herda de `Gauss`.
   - Implementa o método de decomposição LU para decompor a matriz de coeficientes A em matrizes triangulares inferior e superior (L e U), usando eliminação gaussiana.

## Como Usar


#### OBS:
**Para utilizar os métodos é necessário instalar o pacote tabulate com o comando:**   

~~~python
pip install tabulate
~~~

### Importe as classes necessárias:

~~~python
from methods import Jacobi, Seidel, Gauss, LU
~~~

~~~python
# defina a matriz de coeficientes

coefficients = [[a11, a12, a13],
                [a21, a22, a23],
                [a31, a32, a33]]

# defina o vetor de resultados

result = [b1, b2, b3]

# Crie uma instância do Método que deseja utilizar com a matriz de coeficientes e o vetor resultado

method_instance = Method(coefficients, result)
jacobi_instance = Jacobi(coefficients, result)
seidel_instance = Seidel(coefficients, result)
gauss_instance = Gauss(coefficients, result)
lu_instance = LU(coefficients, result)

# Realize a aproximação por um número especificado de iterações e uma tolerância
# Esses parâmetros não são obrigatórios

iterations = 100
tolerance = 1e-6

# Para encontrar as Soluções do sistema por Jacobi ou Gauss-Seidel utilize .approximate()

jacobi_instance.approximate(iterations, tolerance)
seidel_instance.approximate(iterations, tolerance)

# Essas formas também são válidas

jacobi_instance.approximate()
seidel_instance.approximate()

# Para encontrar as soluções por Eliminação de Gauss utilize.gaussian_elimination()

gauss_instance.gaussian_elimination()

# Para encontrar as soluções por Lu utilize .decompose()

lu_instance.decompose()

# É possível recuperar os resultados de todos os métodos para cálculos futuros

jacobi_solutions = jacobi_instance.get_solutions()
seidel_solutions = seidel_instance.get_solutions()
gauss_solutions = gauss_instance.get_solutions()
lu_solutions = lu_instance.get_solutions()
~~~


# Linear System Approximation Methods

This implementation includes five Python classes that offer different methods for approximating solutions to linear systems. These classes are developed by Arthur Andrade D'Olival in October 2023.

## Classes

1. `Method`:
   - Attributes:
     - `coefficients`: Matrix of coefficients for the linear system.
     - `result`: Vector of results for the linear system.

2. `Jacobi`:
   - Inherits from `Method`.
   - Implements the Jacobi method for approximating linear system solutions.

3. `Seidel`:
   - Inherits from `Jacobi`.
   - Implements the Gauss-Seidel method, a variation of the Jacobi method.

4. `Gauss`:
   - Inherits from `Method`
   - Implements the Gaussian elimination method followed by retroactive substitution to find solutions to linear systems.

5. `LU`:
   - Inherits from `Gauss`.
   - Implements the LU decomposition method for breaking down the coefficient matrix A into lower and upper triangular matrices (L and U), using Gaussian elimination.

## How to Use

#### OBS:
**To use the methods it is necessary to install the tabulate package with the following command:**

~~~python
pip install tabulate
~~~

### Import the required classes:
   
~~~python
from methods import Method, Jacobi, Seidel, Gauss, LU
~~~

~~~python
# Define the coefficient matrix

coefficients = [[a11, a12, a13],
                [a21, a22, a23],
                [a31, a32, a33]]

# Define the result vector

result = [b1, b2, b3]

# Create an instance of the Method you want to use with the coefficient matrix and the result vector

method_instance = Method(coefficients, result)
jacobi_instance = Jacobi(coefficients, result)
seidel_instance = Seidel(coefficients, result)
gauss_instance = Gauss(coefficients, result)
lu_instance = LU(coefficients, result)

# Perform the approximation for a specified number of iterations and a tolerance
# These parameters are optional

iterations = 100
tolerance = 1e-6

# To find the solutions of the system using Jacobi or Gauss-Seidel, use .approximate()

jacobi_instance.approximate(iterations, tolerance)
seidel_instance.approximate(iterations, tolerance)

# These forms are also valid

jacobi_instance.approximate()
seidel_instance.approximate()

# To find the solutions using Gaussian Elimination, use .gaussian_elimination()

gauss_instance.gaussian_elimination()

# To find the solutions using LU decomposition, use .decompose()

lu_instance.decompose()

# You can retrieve the results of all methods for future calculations

jacobi_solutions = jacobi_instance.get_solutions()
seidel_solutions = seidel_instance.get_solutions()
gauss_solutions = gauss_instance.get_solutions()
lu_solutions = lu_instance.get_solutions()
~~~
