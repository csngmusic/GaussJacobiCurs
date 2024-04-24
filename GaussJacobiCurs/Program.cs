using System;
using System.Collections.Immutable;
using System.IO;

namespace GaussSeidel_Jacobi
{
    class Program
    {
        static void Main()
        {
            string filePath = "3.txt";
            double[,] coefficients;
            double[] constants;
            // Чтение из файла
            using (StreamReader reader = new StreamReader(filePath))
            {
                int numEquations = int.Parse(reader.ReadLine()); // Чтение количества уравнений
                coefficients = new double[numEquations, numEquations]; // Создание массива коэффициентов
                constants = new double[numEquations]; // Создание массива констант

                // Заполнение массива коэффициентов и массива констант
                for (int i = 0; i < numEquations; i++)
                {
                    string[] coefficients_file = reader.ReadLine().Split(' ');
                    for (int j = 0; j <= numEquations; j++)
                    {
                        if (j != numEquations)
                            coefficients[i, j] = double.Parse(coefficients_file[j]);
                        else
                            constants[i] = double.Parse(coefficients_file[j]);
                    }
                }
            }

            double[] initialGuess = constants; // Начальное приближение
            int iterations = 1000; // Максимальное количество итераций
            double tolerance = 0.0001; // Допустимая погрешность
            
            // Создание экземпляра класса LinearEquationSolver
            LinearEquationSolver solver = new LinearEquationSolver(coefficients, constants, initialGuess, iterations, tolerance);

            // Решение СЛАУ методом Якоби
            Console.WriteLine("Решение методом Якоби:");
            double[] jacobiSolution = solver.JacobiMethod();

            // Решение СЛАУ методом Гаусса-Зейделя
            Console.WriteLine("\nРешение методом Гаусса-Зейделя:");
            double[] gaussSeidelSolution = solver.SolveUsingGaussSeidel();
            
            // Вывод решений СЛАУ
            PrintSolution(jacobiSolution);
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
        private double[,] coefficients; // Матрица коэффициентов
        private double[] constants;     // Вектор констант
        private double[] initialGuess;  // Начальное приближение
        private int iterations;         // Количество итераций
        private double tolerance;       // Допустимая погрешность

        // Конструктор класса
        public LinearEquationSolver(double[,] coefficients, double[] constants, double[] initialGuess, int iterations, double tolerance)
        {
            this.coefficients = coefficients;
            this.constants = constants;
            this.initialGuess = initialGuess;
            this.iterations = iterations;
            this.tolerance = tolerance;
        }

        // Метод решения СЛАУ методом Якоби
        public double[] JacobiMethod()
        {
            int n = constants.Length; // Размер системы уравнений
            double[] currentSolution = new double[n]; // Текущее решение
            double[] nextSolution = new double[n]; // Новое решение
            if (!IsDiagonallyDominant(coefficients))
            {
                RearrangeForDiagonalDominance(ref coefficients, ref constants);
            }
            // Инициализация начальным приближением
            Array.Copy(initialGuess, currentSolution, n);

            for (int iter = 0; iter < iterations; iter++)
            { // Итерации
                for (int i = 0; i < n; i++)
                { // Обновление каждой переменной
                    double sum = constants[i]; // Сумма для i-го уравнения
                    
                    for (int j = 0; j < n; j++)
                    { // Вычитание членов с другими переменными
                        if (j != i)
                        {
                            sum -= coefficients[i, j] * currentSolution[j];
                        }
                    }

                    // Деление на диагональный элемент
                    nextSolution[i] = sum / coefficients[i, i];
                }

                // Проверка сходимости
                if (IsConverged(currentSolution, nextSolution))
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
            if (!IsDiagonallyDominant(coefficients))
            {
                RearrangeForDiagonalDominance(ref coefficients, ref constants);
            }
            int n = constants.Length;
            double[] currentSolution = new double[n];
            double[] nextSolution = new double[n];
            Array.Copy(initialGuess, currentSolution, n);

            for (int iter = 0; iter < iterations; iter++)
            {
                for (int i = 0; i < n; i++)
                {
                    double sum = constants[i];
                    for (int j = 0; j < n; j++)
                    {
                        if (j != i)
                        {
                            sum -= coefficients[i, j] * nextSolution[j];
                        }
                    }
                    nextSolution[i] = sum / coefficients[i, i];
                }

                if (IsConverged(currentSolution, nextSolution))
                {
                    Console.WriteLine($"Метод Гаусса-Зейделя сошелся после {iter + 1} итераций.");
                    return nextSolution;
                }

                Array.Copy(nextSolution, currentSolution, n);
            }

            Console.WriteLine("Метод Гаусса-Зейделя не сошелся за заданное число итераций.");
            return null;
        }
        // Метод для проверки диагональной доминируемости матрицы коэффициентов
        private bool IsDiagonallyDominant(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            // Проверяем каждое уравнение
            for (int i = 0; i < rows; i++)
            {
                double diagonalElement = Math.Abs(matrix[i, i]);
                double sum = 0;

                // Суммируем абсолютные значения всех элементов строки, кроме диагонального элемента
                for (int j = 0; j < cols; j++)
                {
                    if (j != i)
                    {
                        sum += Math.Abs(matrix[i, j]);
                    }
                }

                // Проверяем условие диагональной доминируемости
                if (diagonalElement <= sum)
                {
                    return false; // Матрица не является диагонально доминирующей
                }
            }

            return true; // Матрица является диагонально доминирующей
        }
        private void RearrangeForDiagonalDominance(ref double[,] matrix, ref double[] constants)
        {
            int n = matrix.GetLength(0);
            int maxAttempts = 2 * n; // Максимальное количество попыток перестановок

            // Перебираем каждое уравнение
            for (int i = 0; i < n; i++)
            {
                // Находим сумму абсолютных значений всех элементов строки, кроме диагонального элемента
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += Math.Abs(matrix[i, j]);
                    }
                }

                // Если сумма абсолютных значений меньше или равна диагональному элементу,
                // значит, нужно поменять местами текущее уравнение с другим
                if (sum >= Math.Abs(matrix[i, i]))
                {
                    int attempts = 0; // Счетчик попыток перестановок

                    // Ищем уравнение, с которым текущее уравнение можно поменять местами
                    for (int k = i + 1; k < n && attempts < maxAttempts; k++)
                    {
                        // Меняем местами строки
                        SwapRows(ref matrix, i, k);
                        double temp = constants[i];
                        constants[i] = constants[k];
                        constants[k] = temp;
                        // Проверяем, стала ли матрица диагонально доминирующей после перестановки
                        if (IsDiagonallyDominant(matrix))
                        {
                            return; // Если да, завершаем функцию
                        }

                        // Если нет, отменяем перестановку и пробуем другое уравнение
                        SwapRows(ref matrix, i, k);
                        temp = constants[i];
                        constants[i] = constants[k];
                        constants[k] = temp;
                        attempts++; // Увеличиваем счетчик попыток
                    }
                }
            }
        }

        // Метод для обмена строк матрицы
        private void SwapRows(ref double[,] matrix, int row1, int row2)
        {
            int cols = matrix.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                double temp = matrix[row1, j];
                matrix[row1, j] = matrix[row2, j];
                matrix[row2, j] = temp;
            }
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
