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
   vector<double> grid;          // Исходная сетка
   vector<double> time_grid;     // Исходная сетка по времени

   int nodes_count = 0;          // Общее количество узлов
   int elems_count = 0;          // Общее количество конечных элементов

   vector<double> b;             // Глобальный вектор правой части
   vector<double> solution;      // Полученное решение
   vector<double> true_solution; // Точное решение
   vector<double> xk;            // Вектор начального приближения для МСГ

   // Глобальные матрицы
   Matrix A;                     // Матрица системы
   Matrix G;                     // Матрица жесткости
   Matrix M;                     // Матрица массы

   SLAE slae;                    // СЛАУ
   Test test;                    // Тестовая информация


   HyperbolicalProblem()
   {

   }

   // Функция считывания областей из файла file_name
   // и формирования сетки
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      int n;
      double left, right;

      fin >> left >> right >> n;

      nodes_count = n + 1;
      double h = (right - left) / n;

      grid = vector<double>(nodes_count);
      for(size_t i = 0; i < nodes_count; i++)
         grid[i] = i * h + left;

      elems_count = nodes_count - 1;
   }

   // Функция считывания областей из файла file_name
   // и формирования сетки
   void ReadFormTimeGrid(const string& file_name)
   {
      ifstream fin(file_name);

      int n;
      double left, right;

      fin >> left >> right >> n;

      double h = (right - left) / n;

      time_grid = vector<double>(n + 1);
      for(size_t i = 0; i < n + 1; i++)
         time_grid[i] = i * h + left;
   }
   
   // Инициализация памяти под все необходимые векторы и матрицы
   void InitializeMemory()
   {
      slae = SLAE(nodes_count, 10000, 1e-20);
      b = vector<double>(nodes_count);
      solution = vector<double>(nodes_count);
      true_solution = vector<double>(nodes_count);

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
      // Индекс очередного элемента в треугольнике матрицы
      int to_add_i_tr = 0;
      // Индекс очередного элемента на диагонали матрицы
      int to_add_i_di = 0;

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         double x0 = grid[elem_i];
         double x1 = grid[elem_i + 1];

         double hx = x1 - x0;

         G.diag[to_add_i_di]   += test.lambda / hx * 1.0;
         M.diag[to_add_i_di++] += hx / 6.0 * 2.0;

         G.diag[to_add_i_di] += test.lambda / hx * 1.0;
         M.diag[to_add_i_di] += hx / 6.0 * 2;

         G.bot_tr[to_add_i_tr]   += test.lambda / hx * -1.0;
         M.bot_tr[to_add_i_tr++] += hx / 6.0 * 1.0;
      }

      // Заполнение верхних треуголников
      G.top_tr = G.bot_tr;
      M.top_tr = M.bot_tr;
   }

   // Функция заполнения вектора правой части
   void FillB(const double& t)
   {
      // Индекс очередного элемента на добавление
      int to_add_i = 0;

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         double x0 = grid[elem_i];
         double x1 = grid[elem_i + 1];

         double hx = x1 - x0;

         true_solution[to_add_i] = test.u(x0, t);
         b[to_add_i++]   += hx / 6.0 * (2.0 * test.f(x0, t) + test.f(x1, t));

         true_solution[to_add_i] = test.u(x1, t);
         b[to_add_i] += hx / 6.0 * (test.f(x0, t) + 2 * test.f(x1, t));

      }
   }

   // Функция сборки глобальной матрицы системы
   void AssembleGlobalMatrix()
   {
      for(size_t i = 0; i < nodes_count; i++)
      {
         A.diag[i] = G.diag[i] + test.sigma * M.diag[i] + test.chi * M.diag[i];
      }
      for(size_t i = 0; i < A.tr_size; i++)
      {
         A.bot_tr[i] = G.bot_tr[i] + test.sigma * M.bot_tr[i] + test.chi * M.bot_tr[i];
      }
      A.top_tr = A.bot_tr;
   }

   // Функция учета краевых условий
   void AccountBound(const double& t)
   {
      b[0] = test.u(grid[0], t);
      b[nodes_count - 1] = test.u(grid[nodes_count - 1], t);

      A.diag[0] = 1.0;
      A.diag[nodes_count - 1] = 1.0;

      A.top_tr[0] = 0.0;
      A.bot_tr[A.tr_size - 1] = 0.0;
   }

   // Нахождение решения СЛАУ
   void Solve()
   {
      slae.b = b;
      cout << slae.ConjGradMethod(xk, solution, A) << endl;
   }

   // Вывод решения на временном слое t в поток fout 
   void PrintSolution(ofstream& fout, const double& t)
   {
      fout << "t = " << fixed << t << endl;
      fout << setw(14) << "x";
      fout << setw(14) << "prec" << setw(14) << "calc" << setw(14) << "diff" << setw(5) << "n" << " loc" << endl;

      double norm = 0, norm_u = 0;

      for(int i = 0; i < nodes_count; i++)
      {
         double prec = true_solution[i];
         double calc = solution[i];

         //if(x_i % 32 == 0 && y_i % 32 == 0)
         {
            fout << scientific;
            fout << setw(14) << grid[i];
            fout << setw(14) << prec;
            fout << setw(14) << calc;
            fout << setw(14) << abs(true_solution[i] - solution[i]);
            fout << fixed << setw(5) << i;

            fout << endl;

         }
         norm_u += prec * prec;
         norm += abs(prec - calc) * abs(prec - calc);
      }
   }
};