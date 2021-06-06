#pragma once
#include <vector>
#include <fstream>
using namespace std;

class Matrix
{
public:
   int size = 0;               // Размер матрицы
   int tr_size = 0;            // Количество элементов в треугольнике
   
   vector<int> ind;            // Указатели начала строк
   vector<int> columns_ind;    // Номера столбцов внедиагональных элементов

   vector<double> top_tr;      // Верхний треугольник
   vector<double> bot_tr;      // Нижний треугольник
                               
   vector<double> diag;        // Диагональ

   // Конструктор для матрицы с известным количеством элементов в треуголниках
   Matrix(const int& t_size, const int& t_tr_size) : size(t_size), tr_size(t_tr_size) 
   {
      top_tr = vector<double>(tr_size);
      bot_tr = vector<double>(tr_size);
      columns_ind = vector<int>(tr_size);
      diag = vector<double>(size);
      ind = vector<int>(size + 1);
   }

   Matrix(const Matrix& mat)
   {
      size = mat.size;
      tr_size = mat.tr_size;

      top_tr = mat.top_tr;
      bot_tr = mat.bot_tr;
      diag = mat.diag;
      ind = mat.ind;
      columns_ind = mat.columns_ind;
   }

   Matrix()
   {

   }

   // Функция умножения матрицы на вектор vec, результат в res
   void MatVecMult(const vector<double>& vec, vector<double>& res,
      const vector<double>& bot_tr, const vector<double>& top_tr)
   {
      for (int i = 0; i < size; i++)
         res[i] = 0;

      for (int i = 0; i < size; i++)
      {
         res[i] += vec[i] * diag[i];

         int prof_len = ind[i + 1] - ind[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = ind[i] + k;
            int j = columns_ind[i_in_gg];
            res[i] += vec[j] * bot_tr[i_in_gg];
            res[j] += vec[i] * top_tr[i_in_gg];
         }
      }
   }

   // Функция умножения матрицы на вектор vec, результат в res
   void MatVecMult(const vector<double>& vec, vector<double>& res)
   {
      MatVecMult(vec, res, bot_tr, top_tr);
   }
};