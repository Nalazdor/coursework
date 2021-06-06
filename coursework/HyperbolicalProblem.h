#pragma once
#include <vector>
#include <fstream>
#include <iomanip>
#include "SLAE.h"
#include "Test.h"

using namespace std;

class HyperbolicalProblem
{
public:
   vector<double> grid;             // Сетка по пространству
   vector<double> time_grid;        // Сетка по времени

   int nodes_count = 0;             // Количество узлов
   int elems_count = 0;             // Количество конечных элементов

   vector<double> b;                // Глобальный вектор правой части

   vector<double> solution;         // Полученное решение
   vector<double> true_solution;    // Точное решение

   vector<double> xk;               // Вектор начального приближения для МСГ

   vector<double> prev_solution1;   // Решение на предыдущей итерации по времени
   vector<double> prev_solution2;   // Решение на предпредыдущей итерации по времени

   // Вспомогательные векторы для итерации по времени
   vector<double> Mqj_1, Mqj_2;

   // Глобальные матрицы
   Matrix A;                        // Матрица системы
   Matrix G;                        // Матрица жесткости
   Matrix M;                        // Матрица массы

   SLAE slae;                       // СЛАУ
   Test test;                       // Тестовая информация

   HyperbolicalProblem()
   {

   }

   // Функция считывания и формирования сетки по пространству 
   // из файла file_name
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      int n;              // Количество разбиений
      double left, right; // Левая и правая границы расчетной области

      fin >> left >> right >> n;

      nodes_count = n + 1;
      double h = (right - left) / n; // Шаг разбиения

      grid = vector<double>(nodes_count);
      // Формирование сетки по пространству
      for(size_t i = 0; i < nodes_count; i++)
         grid[i] = i * h + left;

      // Рачет количества конечных элементов
      elems_count = nodes_count - 1;

      fin.close();
   }

   // Функция считывания и формирования сетки по времени 
    // из файла file_name
   void ReadFormTimeGrid(const string& file_name)
   {
      ifstream fin(file_name);

      int n;              // Количество разбиений
      double left, right; // Первая и последняя точки сетки по времени

      fin >> left >> right >> n;

      double h = (right - left) / n; // Шаг разбиения

      time_grid = vector<double>(n + 1);
      // Формирование сетки по времени
      for(size_t i = 0; i < n + 1; i++)
         time_grid[i] = i * h + left;

      fin.close();
   }
   
   // Инициализация памяти под все необходимые векторы и матрицы
   void InitializeMemory()
   {
      slae = SLAE(nodes_count, 10000, 1e-20);
      b = vector<double>(nodes_count);
      solution = vector<double>(nodes_count);
      true_solution = vector<double>(nodes_count);

      prev_solution1 = vector<double>(nodes_count);
      prev_solution2 = vector<double>(nodes_count);

      Mqj_1 = vector<double>(nodes_count);
      Mqj_2 = vector<double>(nodes_count);

      A = Matrix(nodes_count, nodes_count - 1);
      G = Matrix(nodes_count, nodes_count - 1);
      M = Matrix(nodes_count, nodes_count - 1);

      xk = vector<double>(nodes_count);
   }

   // Функция формирования портрета матрицы
   void FormPortrait(Matrix& m)
   {
      m.ind[0] = 0;
      m.ind[1] = 0;

      for(size_t i = 1; i < nodes_count; i++)
      {
         m.ind[i + 1] = m.ind[i] + 1;
         m.columns_ind[i - 1] = i - 1;
      }
   }

   // Функция заполнения матриц жесткости и массы
   void FillMatrices()
   {
      // Индекс очередного элемента на добавление в треугольник матрицы
      int to_add_i_tr = 0;
      // Индекс очередного элемента на добавление в диагональ матрицы
      int to_add_i_di = 0;

      // Цикл по конечным элементам
      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         // Координаты узлов конечного элемента
         double x0 = grid[elem_i];
         double x1 = grid[elem_i + 1];

         // Ширина конечного элемента
         double hx = x1 - x0;

         // Заполнение матрицы жесткости
         G.diag[to_add_i_di]     += test.lambda / hx * 1.0;
         G.diag[to_add_i_di + 1] += test.lambda / hx * 1.0;
         G.bot_tr[to_add_i_tr]   += test.lambda / hx * -1.0;

         // Заполнение матрицы массы
         M.diag[to_add_i_di]     += hx / 6.0 * 2.0;
         M.diag[to_add_i_di + 1] += hx / 6.0 * 2.0;
         M.bot_tr[to_add_i_tr]   += hx / 6.0 * 1.0;

         to_add_i_di++;
         to_add_i_tr++;
      }

      // Заполнение верхних треуголников
      // в силу симметричности матрицы
      G.top_tr = G.bot_tr;
      M.top_tr = M.bot_tr;
   }

   // Функция заполнения вектора правой части
   void FillB(const double& t)
   {
      // Индекс очередного элемента на добавление
      int to_add_i = 0;

      // Цикл по конечным элементам
      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         // Координаты узлов конечного элемента
         double x0 = grid[elem_i];
         double x1 = grid[elem_i + 1];

         // Ширина конечного элемента
         double hx = x1 - x0;

         // Заполнение вектора точного решения
         true_solution[to_add_i]     = test.u(x0, t);
         true_solution[to_add_i + 1] = test.u(x1, t);

         // Заполнение вектора правой части системы
         b[to_add_i]     += hx / 6.0 * (2.0 * test.f(x0, t) + test.f(x1, t));
         b[to_add_i + 1] += hx / 6.0 * (test.f(x0, t) + 2 * test.f(x1, t));

         to_add_i++;
      }
   }

   // Функция сборки глобальной матрицы системы
   void AssembleGlobalMatrix(const double& c1, const double& c2)
   {
      // Сборка диагонали матрицы
      for(size_t i = 0; i < nodes_count; i++)
         A.diag[i] = G.diag[i] + c1 * test.sigma * M.diag[i] + c2 * test.chi * M.diag[i];

      // Сборка нижнего треугольнка матрицы
      for(size_t i = 0; i < A.tr_size; i++)
         A.bot_tr[i] = G.bot_tr[i] + c1 * test.sigma * M.bot_tr[i] + c2 * test.chi * M.bot_tr[i];

      A.top_tr = A.bot_tr;
   }

   // Функция учета краевых условий
   void AccountBound(const double& t)
   {
      // Точное решение в вектор правой части на граничных узлах
      b[0] = test.u(grid[0], t);
      b[nodes_count - 1] = test.u(grid[nodes_count - 1], t);

      // Единица на главной диагонали на граничных узлах
      A.diag[0] = 1.0;
      A.diag[nodes_count - 1] = 1.0;

      // Нули в строках матрицы на граничных узлах
      A.top_tr[0] = 0.0;
      A.bot_tr[A.tr_size - 1] = 0.0;
   }

   // Нахождение решения СЛАУ
   void Solve()
   {
      slae.b = b;
      cout << slae.ConjGradMethod(xk, solution, A) << endl;
   }

   // Базисная функция 1 на шаблонном конечном элементе
   double phi1(const double& xi)
   {
      return 1 - xi;
   }

   // Базисная функция 2 на шаблонном конечном элементе
   double phi2(const double& xi)
   {
      return xi;
   }

   // Вывод решения на временном слое t в поток fout 
   void PrintSolution(ofstream& fout, const double& t)
   {
      // Шапка таблицы
      fout << "t = " << fixed << t << endl;
      fout << setw(14) << "x";
      fout << setw(14) << "prec" << setw(14) << "calc" << setw(14) << "diff" << setw(5) << "n" << endl;

      double norm = 0, norm_u = 0;

      // Цикл по узлам сетки
      for(int i = 0; i < nodes_count; i++)
      {
         double prec = true_solution[i]; // Точное решение в узле
         double calc = solution[i];      // Полученное решение в узле

         //if(i % 32 == 0)
         {
            fout << scientific;
            fout << setw(14) << grid[i];
            fout << setw(14) << prec;
            fout << setw(14) << calc;
            fout << setw(14) << abs(true_solution[i] - solution[i]);
            fout << fixed << setw(5) << i;

            fout << endl;

            // Накопление вектора точного решения
            norm_u += prec * prec;

            // Накопление вектора погрешности
            norm += abs(prec - calc) * abs(prec - calc);
         }
      }

      //// Блок вывода для вывода решения в произвольной точке расчетной области
      //double x = 1.5;
      //double h = grid[1] - grid[0];
      //int elem_i = floor((x - grid[0]) / h);
      //double xi = (x - grid[elem_i]) / h;

      //double calc = phi1(xi) * solution[elem_i] + phi2(xi) * solution[elem_i + 1];
      //double prec = test.u(x, t);

      //fout << scientific;
      //fout << setw(14) << x;
      //fout << setw(14) << prec;
      //fout << setw(14) << calc;
      //fout << setw(14) << abs(prec - calc);
      //fout << fixed << setw(5) << "-";

      //fout << " point" << endl;

      // Расчет и вывод норм векторов относительной и абсолютной погрешности решения
      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*||" << scientific << sqrt(norm) << endl;
   }

   // Итерация трехслойной неявной схемой по времени, вывод в fout
   void IterateTime(ofstream& fout)
   {
      // Сборка матриц, которую можно провести один раз до начала итераций
      FillMatrices();

      // Получение векторов начального приближения на первых двух временных слоях
      for(int x_i = 0; x_i < nodes_count; x_i++)
      {
         prev_solution2[x_i] = test.u(grid[x_i], time_grid[0]);
         prev_solution1[x_i] = test.u(grid[x_i], time_grid[1]);
      }

      // Цикл по временным слоям
      for(int time_i = 2; time_i < time_grid.size(); time_i++)
      {
         // Временные слои
         double t2 = time_grid[time_i - 2];
         double t1 = time_grid[time_i - 1];
         double t0 = time_grid[time_i - 0];

         // Вспомогательные переменные
         double dt  = t0 - t2;
         double dt0 = t0 - t1;
         double dt1 = t1 - t2;

         // Числители дробей
         double den2 = dt1 * dt;
         double den1 = dt1 * dt0;
         double den0 = dt * dt0;

         // Первые произодные
         double n2 = dt0;
         double n1 = dt;
         double n0 = dt + dt0;

         // Вторые производные
         double o2 = 2;
         double o1 = 2;
         double o0 = 2;

         // Обнуление веткора правой части
         for(int i = 0; i < nodes_count; i++)
            b[i] = 0;

         // Заполнение вектора правой части
         FillB(t0);

         M.MatVecMult(prev_solution1, Mqj_1, M.bot_tr, M.top_tr);
         M.MatVecMult(prev_solution2, Mqj_2, M.bot_tr, M.top_tr);

         // Сборка глобальной матрицы системы
         AssembleGlobalMatrix(n0 / den0, o0 / den0);

         // Сборка вектора правой части для трехслойной неявной схемы
         for(int i = 0; i < nodes_count; i++)
            b[i] += - n2 / den2 * test.sigma * Mqj_2[i] + n1 / den1 * test.sigma * Mqj_1[i]
                    - o2 / den2 * test.chi * Mqj_2[i] + o1 / den1 * test.chi * Mqj_1[i];

         // Учет первых краевых условий
         AccountBound(t0);

         // Решение полученной системы
         Solve();

         // Вывод решения
         //if(time_i == time_grid.size() - 1)
         {
            PrintSolution(fout, t0);
            fout << endl;
         }

         // Обновление векторов решения на предыдущих итерациях
         prev_solution2 = prev_solution1;
         prev_solution1 = solution;
      }
   }
};