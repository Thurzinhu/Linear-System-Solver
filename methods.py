from tabulate import tabulate

class Method:
    def __init__(self, coefficients, result):
        self.coefficients = coefficients
        self.result = result
        self._solutions = {}

        for i in range(len(coefficients)):
            if coefficients[i][i] == 0:
                raise ZeroDivisionError("Diagonal coefficients must not be 0")

    @property
    def solutions(self):
        return self._solutions


class Jacobi(Method):
    def __init__(self, coefficients, result):
        super().__init__(coefficients, result)

    
    def get_first_approximation(self) -> list:
        """
            This function determines the first root approximation values
            by dividing the corresponding result value by the diagonal coefficient

            Returns: A list with the first approximated root values
        """
        root_approximation = []
        for i, result in enumerate(self.result):
            root_approximation.append(result / self.coefficients[i][i])

        return root_approximation


    def approximate(self, tolerance=10 ** -2, max_iterations=15):
        """
            Applying Jacobi's algorithm to find the roots.

            tolerance (float, optional): Determines the precision expected. Defaults to 0.01
            max_iterations (int, optional): Determines the maximum number of iterations. Defaults to 15.
        """
        root_approximation = self.get_first_approximation()

        error = []
        for _ in self.coefficients[0]:
            error.append(1)

        print(" Método de Jacobi ".center(70, "#") + '\n')
        # print(self.get_columns())

        table = []
        for interation in range(max_iterations): 
            # print(self.get_values(interation, root_approximation, error))
            new_approximation = []
            for i in range(len(self.coefficients)):
                sum = self.result[i]
                for j in range(len(self.coefficients[i])):
                    if i != j:
                        sum += -self.coefficients[i][j] * root_approximation[j]

                new_approximation.append(sum / self.coefficients[i][i])
                error[i] = abs(new_approximation[i] - root_approximation[i]) / abs(new_approximation[i])

            root_approximation = new_approximation
            table.append(self.get_values(interation, root_approximation, error))
            
            if max(error) < tolerance:
                # print(self.get_values(interation + 1, root_approximation, error))
                break
        
        for i in range(len(self.result)):
            self.solutions[f"x({i})"] = root_approximation[i]
        
        print(tabulate(table, headers=self.get_columns(), tablefmt='rounded_grid'))
        print()


    def get_columns(self):
        """
            function that generates the headers of Jacobi's table
        """
        headers = ["iteração"]
        headers.extend(f"x({i})" for i in range(len(self.coefficients[0])))
        headers.extend(f"err({i})" for i in range(len(self.coefficients[0])))
        headers.append("err")

        return headers


    def get_values(self, interation, values, errors):
        """
            function that returns a list with all Jacobi's table values
        """
        row = [interation]
        row.extend(f"{value:8.6f}" for value in values)
        row.extend(f"{error:8.6f}" for error in errors)
        row.append(f"{max(errors):8.6f}")
        
        return row


class Seidel(Jacobi, Method):
    def __init__(self, coefficients, result):
        super().__init__(coefficients, result)


    def approximate(self, tolerance=10 ** -2, max_iterations=15):
        """
            Applying Seidel's algorithm to find the roots.

            tolerance (float, optional): Determines the precision expected. Defaults to 0.01
            max_iterations (int, optional): Determines the maximum number of iterations. Defaults to 15.
        """
        root_approximation = self.get_first_approximation()

        error = []
        for _ in self.coefficients[0]:
            error.append(1)

        print(" Método de Gauss Seidel ".center(70, "#") + '\n')
        # print(self.get_columns())

        table = []
        for interation in range(max_iterations): 
            # print(self.get_values(interation, root_approximation, error))
            for i in range(len(self.coefficients)):
                sum = self.result[i]
                previous_value = root_approximation[i]
                for j in range(len(self.coefficients[i])):
                    if i != j:
                        sum += -self.coefficients[i][j] * root_approximation[j]
                
                root_approximation[i] = sum / self.coefficients[i][i]
                error[i] = abs(root_approximation[i] - previous_value) / abs(root_approximation[i])

            table.append(self.get_values(interation, root_approximation, error))
            if max(error) < tolerance:
                # print(self.get_values(interation + 1, root_approximation, error))
                break

        for i in range(len(self.result)):
            self.solutions[f"x({i})"] = root_approximation[i]

        print(tabulate(table, headers=self.get_columns(), tablefmt='rounded_grid'))
        print()


class Gauss(Method):
    def __init__(self, coefficients, result):
        super().__init__(coefficients, result)

        # keeping original values to calculate the residual after defining the roots
        self.original_coefficients = [row[:] for row in coefficients]
        self.original_result = result[:]

        identity = []
        for i in range(len(self.coefficients)):
            row = []
            for j in range(len(self.result)):
                row.append(1) if i == j else row.append(0)
            identity.append(row)
        
        # defining permutation and lower triangular matrices
        self.permutations = [row[:] for row in identity]
        self.lower_triangular_matrix = [row[:] for row in identity]


    def gaussian_elimination(self):
        """
            Applying gaussian elimination by taking every value below the diagonal 
            and adding with a multiple of the pivot line
        """
        print(" Eliminação de Gauss ".center(70, "#") + '\n')
        self.print_values()

        for i in range(len(self.result) - 1):
            modified_rows = {}
            self.find_row_with_largest_coefficient(i) # chaging rows based on largest coefficient
            for j in range(i + 1, len(self.coefficients)):
                if self.coefficients[j][i] != 0:
                    pivot = -(self.coefficients[j][i] / self.coefficients[i][i])
                    modified_rows[f"L({j})"] = pivot
                    self.lower_triangular_matrix[j][i] = pivot

                    for k in range(len(self.coefficients[j])):
                        self.coefficients[j][k] += (pivot * self.coefficients[i][k])
                    self.result[j] += (pivot * self.result[i])
            
            self.print_values(i, modified_rows)


    def print_values(self, interation=-1, modified_rows={}):
        if interation == -1:
            print("Matriz original".center(70, '-') + '\n')
        else:
            for key, value in modified_rows.items():
                print(f"{key} <- {key} + L({interation}) * {value:4.4f}")
            print()

            print(f"Matriz iteração: {interation + 1}".center(70, '-') + '\n')

        for i in range (len(self.coefficients)):
            print("[", end='')
            for j in range(len(self.result)):
                print(f"{self.coefficients[i][j]:8.4f}", end='')
            print(f" ] [ {self.result[i]:8.4f} ]")
        
        print()


    def retroactive_substitution(self):
        """
            Creating a dictionary containing variable associated with its value
            Applying substitution from bottom-up
        """
        for i in range(len(self.coefficients) - 1, -1, -1):
            sum = 0
            for j in range(i + 1, len(self.coefficients)):
                sum += self.solutions[f"x({j})"] * self.coefficients[i][j]
            
            self.solutions[f"x({i})"] = (self.result[i] - sum) / self.coefficients[i][i]

        table = [[key, value] for key, value in self.solutions.items()]

        print(" Raízes de X ".center(70, "#") + '\n')
        print(tabulate(table, headers=["Variável", "Valor"], tablefmt="rounded_grid"))
        print()
        # for i in range(len(self.result)):
            # print(f"x({i}) = {roots[f'x({i})']:8.4f}")


        self.get_remainder()


    def get_remainder(self):
        """
            Calculating residual by substituting the roots on each line of original matrix
            and comparing with the expected result
        """
        table = []
        for i in range(len(self.original_coefficients)):
            sum = 0
            for j in range(len(self.original_result)):
                sum += self.original_coefficients[i][j] * self.solutions[f"x({j})"]

            table.append([f"{self.original_result[i]:8.4f}", f"{sum:8.4f}", f"{self.original_result[i] - sum:8.4f}"])
            # table.append(sum)

        print(" Resíduo ".center(70, "#") + '\n')
        headers = ["Val. Esperado", "Val. obtido", "Resíduo"]
        print(tabulate(table, headers=headers, tablefmt='rounded_grid'))
        print()

        # print("Val. Esperado    Val. obtido    Resíduo")
        # for i in range(len(table)):
        #     print(f"| {self.original_result[i]:13.4f} |", end='')
        #     print(f"{table[i]:11.4f} |", end='')
        #     print(f"{self.original_result[i] - table[i]:8.4f} |")
        # print()


    def find_row_with_largest_coefficient(self, iteration):
        """
            Changing rows based on largest absolute coefficient value
            Updating coefficients, result an permutations matrices
        """
        column = iteration
        largest_coefficient = abs(self.coefficients[iteration][iteration])
        largest_coefficient_row = iteration

        for i in range(iteration + 1, len(self.result)):
            if abs(self.coefficients[i][column]) > largest_coefficient:
                largest_coefficient_row = i
                largest_coefficient = self.coefficients[i][column]

        if largest_coefficient_row != iteration:
            self.coefficients[iteration], self.coefficients[largest_coefficient_row] = (
                self.coefficients[largest_coefficient_row],
                self.coefficients[iteration]
            )
            self.permutations[iteration], self.permutations[largest_coefficient_row] = (
                self.permutations[largest_coefficient_row],
                self.permutations[iteration]
            )
            self.result[iteration], self.result[largest_coefficient_row] = (
                self.result[largest_coefficient_row], 
                self.result[iteration]
            )

            print('\n' + f"L({iteration}) <- L({largest_coefficient_row})")
            print(f"L({largest_coefficient_row}) <- L({iteration})" + '\n')


class LU(Gauss, Method):
    def __init__(self, coefficients, result):
        super().__init__(coefficients, result)


    def decompose(self):
        """
            Applying gaussian elimination and getting lower triangular, permutations,
            B * P and Y matrices
        
            A * X = P * B
            (L * U) * X = P * B
        
            (1) U * X = Y
            (2) L * Y = P * B
        """
        super().gaussian_elimination()

        print(" Método da decomposição LU ".center(70, "#") + '\n')

        self.permuted_result = self.get_permuted_resul()

        self.print_matrices()
        self.result = self.get_Y()


    def get_Y(self):
        """
            Calculating Y matrix
            L * Y = P * B

            Returns: Y 
        """
        roots = {}

        for i in range(len(self.coefficients)):
            sum = 0
            for j in range(i - 1, -1, -1):
                sum += self.lower_triangular_matrix[i][j] * roots[f"y({j})"]
            roots[f"y({i})"] = (self.permuted_result[i] - sum) / self.lower_triangular_matrix[i][i]

            table = [[key, value] for key, value in roots.items()]

        print(" Raízes de Y ".center(70, "#") + '\n')
        print(tabulate(table, headers=["Variável", "Valor"], tablefmt="rounded_grid"))
        print()

        return [value for value in roots.values()]


    def get_permuted_resul(self):
        """
            Multiplying B by P

            returns: permuted result 
        """
        permuted_result = []

        for i in range(len(self.permutations)):
            sum = 0
            for j in range(len(self.original_result)):
                sum += self.permutations[i][j] * self.original_result[j]
            permuted_result.append(sum)

        return permuted_result


    def print_matrices(self):
        matrices = {
            " Matriz Triangular Inferior(L)" : self.lower_triangular_matrix, 
            " Matriz Triangular Superior(U)" : self.coefficients, 
            " Matriz das Permutações" : self.permutations,
        }

        for name, matrix in matrices.items():
            print(name.center(70, '-') + '\n')
            for i in range(len(self.coefficients)):
                print("[", end='')
                for j in range(len(self.result)):
                    print(f"{matrix[i][j]:8.4f}", end=' ')
                print("]")
            print()

        print(" Matriz P * B ".center(70, '-') + '\n')
        for i in range(len(self.result)):
            print(f"[ {self.permuted_result[i]:8.4f} ]")
        print()