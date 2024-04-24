using System;
using System.Collections.Immutable;
using System.IO;

namespace GaussSeidel_Jacobi
{
    class Program
    {
        static void Main()
        {
            string filePath = "1.txt";
            double[,] A;
            double[] b;
            // Чтение из файла
            using (StreamReader reader = new StreamReader(filePath))
            {
                int numEquations = int.Parse(reader.ReadLine()); // Чтение количества уравнений
                A = new double[numEquations, numEquations]; // Создание массива коэффициентов
                b = new double[numEquations]; // Создание массива констант
                // Заполнение массива коэффициентов и массива констант
                for (int i = 0; i < numEquations; i++)
                {
                    string[] coefficients_file = reader.ReadLine().Split(' ');
                    for (int j = 0; j <= numEquations; j++)
                    {
                        if (j != numEquations)
                            A[i, j] = double.Parse(coefficients_file[j]);
                        else
                            b[i] = double.Parse(coefficients_file[j]);
                    }
                }
            }
            double[,] B = new double[b.Length, b.Length];
            double[,] E = new double[b.Length, b.Length];
            double[,] D = new double[b.Length, b.Length];
            double[,] L = new double[b.Length, b.Length];
            double[] g = new double[b.Length];
            for (int i = 0; i < b.Length; i++)
            {
                for (int j = 0; j < b.Length; j++)
                {
                    E[i, j] = 0;
                    D[i, j] = 0;
                    if (i == j)
                    {
                        E[i, j] = 1;
                        D[i, j] = A[i, j];
                    }
                    if (i > j)
                        L[i, j] = A[i, j];
                }
                g[i] = 0;
            }
            
            double[] initialGuess = b; // Начальное приближение
            int iterations = 10000; // Максимальное количество итераций
            double tolerance = 0.00001; // Допустимая погрешность

            // Создание экземпляра класса LinearEquationSolver
            LinearEquationSolver solver = new LinearEquationSolver(A, B, E, D, L, b, g, initialGuess, iterations, tolerance);

            // Решение СЛАУ методом Якоби
            Console.WriteLine("Решение методом Якоби:");
            double[] jacobiSolution = solver.JacobiMethod();
            PrintSolution(jacobiSolution);

            // Решение СЛАУ методом Гаусса-Зейделя
            Console.WriteLine("\nРешение методом Гаусса-Зейделя:");
            double[] gaussSeidelSolution = solver.SolveUsingGaussSeidel();
            PrintSolution(gaussSeidelSolution);

            Console.ReadLine();
        }
        // Метод для вывода решения на экран
        static void PrintSolution(double[] solution)
        {
            if (solution != null)
            {
                Console.Write("Решение: ");
                foreach (var value in solution)
                {
                    Console.Write($"{value:F4} ");
                }
                Console.WriteLine();
            }
            else
            {
                Console.WriteLine("Решение не найдено.");
            }
        }
    }
    class LinearEquationSolver
    {
        private double[,] A; // Матрица коэффициентов
        private double[] b;     // Вектор констант
        private double[,] B;
        private double[,] E;
        private double[,] D;
        private double[,] L;
        private double[] g;
        private double[] initialGuess;  // Начальное приближение
        private int iterations;         // Количество итераций
        private double tolerance;       // Допустимая погрешность

        // Конструктор класса
        public LinearEquationSolver(double[,] A, double[,] B, double[,] E, double[,] D, double[,] L, double[] b, double[] g,
                                    double[] initialGuess, int iterations, double tolerance)
        {
            this.A = A;
            this.B = B;
            this.E = E;
            this.D = D;
            this.L = L;
            this.b = b;
            this.g = g;
            this.initialGuess = initialGuess;
            this.iterations = iterations;
            this.tolerance = tolerance;
        }

        // Метод решения СЛАУ методом Якоби
        public double[] JacobiMethod()
        {
            bool check = false;
            int n = b.Length; // Размер системы уравнений
            double[] currentSolution = new double[n]; // Текущее решение
            double[] nextSolution = new double[n]; // Новое решение
            double[,] D_1 = new double[n, n];
            // Инициализация начальным приближением
            Array.Copy(initialGuess, currentSolution, n);
            Array.Copy(D, D_1, n * n);
            D_1 = InverseMatrix(D_1);
            for (int i = 0; i < n; i++)
                g[i] = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B[i, j] = E[i, j];
                    for (int k = 0; k < n; k++)
                        B[i, j] -= D_1[i, k] * A[k, j];
                    g[i] += D_1[i, j] * b[j];
                }
            }
            double norm = CalculateMatrix1Norm(B);
            if (norm > 1)
            {
                Console.WriteLine("Норма матрицы равна: " + norm);
                Console.WriteLine("Метод может не сойтись!");
            }
            else
                Console.WriteLine("Норма матрицы равна: " + norm);
            for (int iter = 0; iter < iterations; iter++)
            { // Итерации
                for (int i = 0; i < n; i++)
                { // Обновление каждой переменной
                    nextSolution[i] = g[i];
                    for (int j = 0; j < n; j++)
                        nextSolution[i] += B[i, j] * currentSolution[j];
                }
                foreach(var solution in nextSolution)
                {
                    if(solution == Double.NegativeInfinity || Double.IsNaN(solution))
                    {
                        check = true;
                        Console.WriteLine($"Метод Якоби не сошелся после {iter + 1} итераций.");
                        return null;
                    }
                }    
                // Проверка сходимости
                if (IsConverged(currentSolution, nextSolution) && !check)
                {
                    Console.WriteLine($"Метод Якоби сошелся после {iter + 1} итераций.");
                    return nextSolution;
                }

                // Переход к следующему шагу итерации
                Array.Copy(nextSolution, currentSolution, n);
            }

            Console.WriteLine("Метод Якоби не сошелся за заданное число итераций.");
            return null;
        }
        // Метод решения СЛАУ методом Гаусса-Зейделя
        public double[] SolveUsingGaussSeidel()
        {
            bool check = false;
            int n = b.Length; // Размер системы уравнений
            double[] currentSolution = new double[n]; // Текущее решение
            double[] nextSolution = new double[n]; // Новое решение
            double[,] D_2 = new double[n, n];
            // Инициализация начальным приближением
            Array.Copy(initialGuess, currentSolution, n);
            Array.Copy(D, D_2, n*n);
            for (int i = 0; i < n; i++)
            {
                g[i] = 0;
                for(int j = 0; j < n; j++)
                {
                    D_2[i, j] += L[i, j];
                }
            }
            D_2 = InverseMatrix(D_2);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B[i, j] = E[i, j];
                    for (int k = 0; k < n; k++)
                        B[i, j] -= D_2[i, k] * A[k, j];
                    g[i] += D_2[i, j] * b[j];
                }
            }
            double norm = CalculateMatrix1Norm(B);
            if (norm > 1)
            {
                Console.WriteLine("Норма матрицы равна: " + norm);
                Console.WriteLine("Метод может не сойтись!");
            }
            else
                Console.WriteLine("Норма матрицы равна: " + norm);
            for (int iter = 0; iter < iterations; iter++)
            { // Итерации
                for (int i = 0; i < n; i++)
                { // Обновление каждой переменной
                    nextSolution[i] = g[i];
                    for (int j = 0; j < n; j++)
                        nextSolution[i] += B[i, j] * currentSolution[j];
                }
                foreach (var solution in nextSolution)
                {
                    if (solution == Double.NegativeInfinity || Double.IsNaN(solution))
                    {
                        check = true;
                        Console.WriteLine($"Метод Гаусса-Зейделя не сошелся после {iter + 1} итераций.");
                        return null;
                    }
                }
                // Проверка сходимости
                if (IsConverged(currentSolution, nextSolution) && !check)
                {
                    
                    Console.WriteLine($"Метод Гаусса-Зейделя сошелся после {iter + 1} итераций.");
                    return nextSolution;
                }

                // Переход к следующему шагу итерации
                Array.Copy(nextSolution, currentSolution, n);
            }

            Console.WriteLine("Метод Гаусса-Зейделя не сошелся за заданное число итераций.");
            return null;
        }
        private double[,] InverseMatrix(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] augmentedMatrix = AugmentMatrix(matrix);

            // Applying Gaussian elimination
            for (int i = 0; i < n; i++)
            {
                // Finding pivot row
                int pivotRow = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(augmentedMatrix[j, i]) > Math.Abs(augmentedMatrix[pivotRow, i]))
                    {
                        pivotRow = j;
                    }
                }

                // Swapping rows if necessary
                if (pivotRow != i)
                {
                    for (int k = 0; k < 2 * n; k++)
                    {
                        double temp = augmentedMatrix[i, k];
                        augmentedMatrix[i, k] = augmentedMatrix[pivotRow, k];
                        augmentedMatrix[pivotRow, k] = temp;
                    }
                }
                // Normalizing pivot row
                double pivotValue = augmentedMatrix[i, i];
                for (int k = 0; k < 2 * n; k++)
                {
                    augmentedMatrix[i, k] /= pivotValue;
                }
                // Elimination
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        double factor = augmentedMatrix[j, i];
                        for (int k = 0; k < 2 * n; k++)
                        {
                            augmentedMatrix[j, k] -= factor * augmentedMatrix[i, k];
                        }
                    }
                }
            }
            // Extracting inverse matrix
            double[,] inverse = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverse[i, j] = augmentedMatrix[i, j + n];
                }
            }

            return inverse;
        }
        private double[,] AugmentMatrix(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] augmentedMatrix = new double[n, 2 * n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    augmentedMatrix[i, j] = matrix[i, j];
                }
                augmentedMatrix[i, i + n] = 1;
            }
            return augmentedMatrix;
        }
        private double CalculateMatrix1Norm(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            double maxSum = 0;

            // Суммируем абсолютные значения элементов по столбцам и выбираем максимальную сумму
            for (int j = 0; j < cols; j++)
            {
                double columnSum = 0;
                for (int i = 0; i < rows; i++)
                {
                    columnSum += Math.Abs(matrix[i, j]);
                }
                maxSum = Math.Max(maxSum, columnSum);
            }

            return maxSum;
        }
        // Метод проверки сходимости
        private bool IsConverged(double[] previousSolution, double[] currentSolution)
        {
            int n = previousSolution.Length;
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(currentSolution[i] - previousSolution[i]) > tolerance)
                {
                    return false;
                }
            }
            return true;
        }
    }

}
