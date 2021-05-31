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

   vector<double> solution;      // Полученное решение
   vector<double> true_solution; // Точное решение

   Matrix G, M;                  // Матрицы жесткости и массы
   SLAE slae;                    // СЛАУ
   Test test;                    // Тестовая информация


   // Вспомогательная матрица для построения
   // матрицы массы конечного элемента
   vector<vector<int>> C;

   HyperbolicalProblem()
   {
      C = {
         {4, 2, -1},
         {2, 16, 2},
         {-1, 2, 4}
      };
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

      elems_count = nodes_count / 2;
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

   //// Функция инициализации матрицы
   //void InitMatrix(Matrix& m)
   //{
   //   m.size = 0;
   //   m.bot_tr = vector<double>(elems_count * 3);
   //   m.top_tr = vector<double>(elems_count * 3);

   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      m.size += regions[reg_i].elems_count * 2 + 1;

   //      if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] == grid[(reg_i + 1) * 2])
   //         m.size--;
   //   }

   //   m.di = vector<double>(m.size);
   //   m.ind = vector<int>(m.size + 1);
   //}

   //// Функция формирования портрета глобальной матрицы
   //void FormPortrait(Matrix& m)
   //{
   //   m.ind[0] = 0;
   //   m.ind[1] = 0;
   //   m.ind[2] = 1;
   //   m.ind[slae.m.size] = m.top_tr.size();

   //   int global_i = 3;
   //   int val = 3;

   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      Region* r = &regions[reg_i];

   //      for(int elem_i = reg_i == 0 ? 1 : 0; elem_i < r->elems_count; elem_i++, global_i += 2, val += 2)
   //      {
   //         m.ind[global_i] = val++;
   //         m.ind[global_i + 1] = val;
   //      }

   //      if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] != grid[(reg_i + 1) * 2])
   //         m.ind[global_i++] = val;
   //   }
   //}

   //// Функция заполнения матриц жесткости и массы
   //void FillMatrices(const double& t)
   //{
   //   // Индекс очередного элемента в треугольнике матрицы
   //   int to_add_i_tr = 0;
   //   // Индекс очередного элемента на диагонали матрицы
   //   int to_add_i_di = 0;

   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      Region* r = &regions[reg_i];

   //      for(int elem_i = 0; elem_i < r->elems_count; elem_i++)
   //      {
   //         // Индекс первого узла элемента
   //         int elem_beg_i = elem_i * 2;

   //         // Координаты узлов
   //         vector<double> x_elem = { r->nodes[elem_beg_i], r->nodes[elem_beg_i + 1], r->nodes[elem_beg_i + 2] };
   //         double h = x_elem[2] - x_elem[0];

   //         // Индекс первого узла элемента в глобальной нумерации
   //         int elem_beg_i_glob = elem_i * 2 + r->first_i;

   //         vector<double> q_elem = { u[elem_beg_i_glob], u[elem_beg_i_glob + 1], u[elem_beg_i_glob + 2] };

   //         vector<double> lambda = { test.lambda(q_elem[0], x_elem[0], t), test.lambda(q_elem[1], x_elem[1], t), test.lambda(q_elem[2], x_elem[2], t) };

   //         // Заполнение диагонали матрицы жесткости
   //         G.di[to_add_i_di++] += (lambda[0] * 37 / 30 + lambda[1] * 1.2 - lambda[2] * 0.1) / h;
   //         G.di[to_add_i_di++] += (lambda[0] * 1.6 + lambda[1] * 32 / 15 + lambda[2] * 1.6) / h;
   //         G.di[to_add_i_di]   += (-lambda[0] * 0.1 + lambda[1] * 1.2 + lambda[2] * 37 / 30) / h;

   //         // Заполнение нижнего треугольника матрицы жесткости
   //         G.bot_tr[to_add_i_tr++] += (-lambda[0] * 22 / 15 - lambda[1] * 16 / 15 - lambda[2] * 2 / 15) / h;
   //         G.bot_tr[to_add_i_tr++] += (lambda[0] * 7 / 30 - lambda[1] * 2 / 15 + lambda[2] * 7 / 30) / h;
   //         G.bot_tr[to_add_i_tr++] += (-lambda[0] * 2 / 15 - lambda[1] * 16 / 15 - lambda[2] * 22 / 15) / h;

   //         to_add_i_di -= 2;
   //         to_add_i_tr -= 3;

   //         // Заполнение диагонали матрицы массы
   //         M.di[to_add_i_di++] += test.sigma() * h / 30.0 * C[0][0];
   //         M.di[to_add_i_di++] += test.sigma() * h / 30.0 * C[1][1];
   //         M.di[to_add_i_di]   += test.sigma() * h / 30.0 * C[2][2];

   //         // Заполнение нижнего треугольника массы
   //         M.bot_tr[to_add_i_tr++] += test.sigma() * h / 30.0 * C[1][0];
   //         M.bot_tr[to_add_i_tr++] += test.sigma() * h / 30.0 * C[2][0];
   //         M.bot_tr[to_add_i_tr++] += test.sigma() * h / 30.0 * C[2][1];

   //         vector<double> f_elem = { test.f(x_elem[0], t), test.f(x_elem[1], t), test.f(x_elem[2], t) };

   //         // Заполнение вектора правой части
   //         slae.b[to_add_i_di - 2] += h / 30.0 * (C[0][0] * f_elem[0] + C[0][1] * f_elem[1] + C[0][2] * f_elem[2]);
   //         slae.b[to_add_i_di - 1] += h / 30.0 * (C[1][0] * f_elem[0] + C[1][1] * f_elem[1] + C[1][2] * f_elem[2]);
   //         slae.b[to_add_i_di - 0] += h / 30.0 * (C[2][0] * f_elem[0] + C[2][1] * f_elem[1] + C[2][2] * f_elem[2]);
   //      }
   //   }

   //   // Заполнение верхних треуголников
   //   G.top_tr = G.bot_tr;
   //   M.top_tr = M.bot_tr;
   //}
   //
   //// Функция сбора глобальной матрицы системы
   //void GlobalMatrixAssemble()
   //{
   //   int tr_size = slae.m.bot_tr.size();
   //   for(int i = 0; i < tr_size; i++)
   //      slae.m.bot_tr[i] = G.bot_tr[i] + M.bot_tr[i];

   //   slae.m.top_tr = slae.m.bot_tr;

   //   for(int i = 0; i < slae.m.size; i++)
   //      slae.m.di[i] = G.di[i] + M.di[i];
   //}

   //// Функция учета краевых условий
   //void AccountBound(const double& t)
   //{
   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      Region* r = &regions[reg_i];

   //      if(r->left_bord == 1)
   //      {
   //         /*slae.m.di[r->first_i] = big_num;
   //         slae.b[r->first_i] = big_num * test.u_prec(r->nodes[0], t);*/

   //         slae.m.di[r->first_i] = 1.0;
   //         slae.b[r->first_i] = test.u_prec(r->nodes[0], t);

   //         slae.m.top_tr[slae.m.ind[r->first_i + 1]] = 0;
   //         slae.m.top_tr[slae.m.ind[r->first_i + 2]] = 0;
   //      }

   //      if(r->right_bord == 1)
   //      {
   //         /*slae.m.di[r->first_i + r->nodes_count - 1] = big_num;
   //         slae.b[r->first_i + r->nodes_count - 1] = big_num * test.u_prec(r->nodes[r->nodes_count - 1], t);*/

   //         slae.m.di[r->first_i + r->nodes_count - 1] = 1.0;
   //         slae.b[r->first_i + r->nodes_count - 1] = test.u_prec(r->nodes[r->nodes_count - 1], t);

   //         slae.m.bot_tr[slae.m.ind[r->first_i + r->nodes_count - 1]] = 0;
   //         slae.m.bot_tr[slae.m.ind[r->first_i + r->nodes_count - 1] + 1] = 0;
   //      }
   //   }
   //}

   //// Функция расчета нормы вектора
   //static double CalcNorm(const vector<double>& vec)
   //{
   //   double res = 0;
   //   int size = vec.size();

   //   for(int i = 0; i < size; i++)
   //      res += vec[i] * vec[i];
   //   
   //   return sqrt(res);
   //}

   //// Линеаризация методом простой итераций, вывод итераций в файл file_name
   //void SimpleIterations(const double& t, const double& eps, const double& delta, const int& max_iter, ofstream& fout)
   //{
   //   int n_iter = 0;
   //   double eps_residual = 1;
   //   double delta_residual = 1;

   //   do
   //   {
   //      // Инициализация матриц
   //      InitMatrix(slae.m);
   //      InitMatrix(G);
   //      InitMatrix(M);
   //      slae.b = vector<double>(slae.m.size);

   //      // Сборка матриц жесткости и массы
   //      FormPortrait(slae.m);
   //      FormPortrait(G);
   //      FormPortrait(M);
   //      FillMatrices(t);

   //      // Сборка матрицы системы
   //      GlobalMatrixAssemble();
   //      AccountBound(t);

   //      vec_1 = vector<double>(nodes_count);
   //      slae.m.MatVecMult(u, vec_1);

   //      for(int i = 0; i < nodes_count; i++)
   //         vec_1[i] -= slae.b[i];

   //      eps_residual = CalcNorm(vec_1) / CalcNorm(slae.b);

   //      if(eps_residual < eps)
   //      {
   //         break;
   //      }
   //      else
   //      {
   //         slae.LUDecomp();
   //         slae.ForwardSolver();
   //         slae.BackwardSolver();

   //         for(int i = 0; i < nodes_count; i++)
   //            vec_1[i] = slae.b[i] - u[i];

   //         delta_residual = CalcNorm(vec_1) / CalcNorm(slae.b);

   //         if(delta_residual < delta)
   //         {
   //            break;
   //         }

   //         u = slae.b;

   //         fout << endl << "Non-linear iteration " << n_iter << endl;
   //         PrintSolution(u, t, fout);
   //         n_iter++;
   //      }
   //   } while(0 < n_iter && n_iter <  max_iter);
   //}

   //// Неявная разностная схема по времени, вывод итераций в поток fout
   //void ExplicitScheme(const double& eps, const double& delta, const int& max_iter, const string& file_name)
   //{
   //   ofstream fout(file_name);

   //   // Расчет вектора приближения при t=0
   //   q_1.resize(nodes_count);

   //   int to_add_i = 0;
   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      Region* r = &regions[reg_i];

   //      for(int elem_i = 0; elem_i < r->elems_count; elem_i++)
   //      {
   //         q_1[to_add_i++] = test.u_prec(r->nodes[elem_i * 2], 0.0);
   //         q_1[to_add_i++] = test.u_prec(r->nodes[elem_i * 2 + 1], 0.0);
   //         q_1[to_add_i] = test.u_prec(r->nodes[elem_i * 2 + 2], 0.0);
   //      }
   //   }

   //   for(int time_i = 1; time_i < time_grid.size(); time_i++)
   //   {
   //      for(int i = 0; i < 100; i++)
   //         fout << "#";
   //      fout << endl;

   //      const double t = time_grid[time_i];
   //      const double dt = t - time_grid[time_i - 1];

   //      fout << "t = " << fixed << t << endl;

   //      int n_iter = 0;
   //      double eps_residual = 1;
   //      double delta_residual = 1;

   //      do
   //      {
   //         // Инициализация матриц
   //         InitMatrix(slae.m);
   //         InitMatrix(G);
   //         InitMatrix(M);
   //         slae.b = vector<double>(slae.m.size);
   //         vec_2 = vector<double>(slae.m.size);

   //         // Сборка матриц жесткости и массы
   //         FormPortrait(slae.m);
   //         FormPortrait(G);
   //         FormPortrait(M);
   //         FillMatrices(t);

   //         // Сборка матрицы системы
   //         const int tr_size = slae.m.bot_tr.size();
   //         for(int i = 0; i < tr_size; i++)
   //            slae.m.bot_tr[i] = G.bot_tr[i] + M.bot_tr[i] / dt;

   //         slae.m.top_tr = slae.m.bot_tr;

   //         for(int i = 0; i < slae.m.size; i++)
   //            slae.m.di[i] = G.di[i] + M.di[i] / dt;

   //         // Ðàñ÷åò âåêòîðà ïðàâîé ÷àñòè
   //         M.MatVecMult(q_1, vec_2);

   //         for(int i = 0; i < slae.m.size; i++)
   //            slae.b[i] += vec_2[i] / dt;

   //         AccountBound(t);

   //         vec_1 = vector<double>(nodes_count);
   //         slae.m.MatVecMult(u, vec_1);

   //         for(int i = 0; i < nodes_count; i++)
   //            vec_1[i] -= slae.b[i];

   //         eps_residual = CalcNorm(vec_1) / CalcNorm(slae.b);

   //         if(eps_residual < eps)
   //         {
   //            break;
   //         }
   //         else
   //         {
   //            slae.LUDecomp();
   //            slae.ForwardSolver();
   //            slae.BackwardSolver();

   //            for(int i = 0; i < nodes_count; i++)
   //               vec_1[i] = slae.b[i] - u[i];

   //            delta_residual = CalcNorm(vec_1) / CalcNorm(slae.b);

   //            if(delta_residual < delta)
   //            {
   //               break;
   //            }

   //            u = slae.b;

   //            //fout << endl << "Non-linear iteration " << n_iter << endl;
   //            //PrintSolution(u, t, fout);
   //            n_iter++;
   //         }
   //      } while(0 < n_iter && n_iter < max_iter);

   //      q = u;
   //      q_1 = q;
   //      //if (time_i % 8 == 0) 
   //      {
   //         fout << "Number of iterations: " << n_iter << endl;
   //         fout << "Final result on timestep:" << endl;
   //         PrintSolution(q, t, fout);
   //      }
   //   }

   //   fout.close();
   //}

   //// Вывод решения на основе вектора текущего приближения solution
   //// во временной точке t в поток fout
   //void PrintSolution(const vector<double>& solution, const double& t, ofstream& fout)
   //{
   //   double norm = 0., norm_u = 0.;

   //   fout << "   x              calc           prec      dif            N  location" << endl << fixed;

   //   for(int reg_i = 0; reg_i < regions_count; reg_i++)
   //   {
   //      Region* r = &regions[reg_i];

   //      for(int node_i = 0; node_i < r->nodes_count; node_i++)
   //      {
   //         //if (node_i % 16 == 0)
   //         {
   //            // Индекс узла в глобальной нумерации
   //            int elem_beg_i = node_i + r->first_i;

   //            fout << setw(11) << r->nodes[node_i];
   //            double help_1 = solution[elem_beg_i];
   //            fout << setw(15) << help_1;
   //            double help_2 = test.u_prec(r->nodes[node_i], t);
   //            fout << setw(15) << help_2;
   //            fout << setw(14) << scientific << abs(help_1 - help_2) << fixed;

   //            fout << setw(4) << elem_beg_i << " ";


   //            if (node_i == 0 || node_i == r->nodes_count - 1)
   //               fout << "  border";
   //            else
   //               fout << "  inner";

   //            fout << endl;

   //            norm_u += help_2 * help_2;
   //            norm += abs(help_1 - help_2) * abs(help_1 - help_2);
   //         }
   //      }
   //      PrintSolutionInside(&regions[reg_i], solution, t, fout, norm_u, norm);
   //   }
   //   fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
   //   fout << "||u-u*|| = " << scientific << sqrt(norm) << endl << endl;
   //}

   //// Вывод решения на основе вектора текущего приближения solution
   //// во временной точке t в поток fout
   //void PrintSolutionInside(Region *r, const vector<double>& solution, const double& t, ofstream& fout, double &norm_u, double &norm)
   //{
   //   // Индекс узла в глобальной нумерации
   //   int i = r->first_i;
   //   double h = (r->nodes[2] - r->nodes[0]);
   //   double x = r->nodes[0] + h / 4;

   //   double e = (x - r->nodes[0]) / h;
   //   double fi1 = 2 * (e - 0.5) * (e - 1);
   //   double fi2 = - 4 * e * (e - 1);
   //   double fi3 = 2 * (e - 0.5) * e;
   //   double res = fi1 * solution[i] + fi2 * solution[i + 1] + fi3 * solution[i + 2];

   //   fout << setw(11) << x;
   //   fout << setw(15) << res;
   //   double help_2 = test.u_prec(x, t);
   //   fout << setw(15) << help_2;
   //   fout << setw(14) << scientific << abs(res - help_2) << fixed;
   //   fout << setw(5) << "- ";
   //   fout << "  point" << endl;

   //   norm_u += help_2 * help_2;
   //   norm += abs(res - help_2) * abs(res - help_2);
   //}
};