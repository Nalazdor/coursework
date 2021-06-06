#pragma once
#include "Vector.h"
#include "Matrix.h"
using namespace std;

class SLAE
{
public:
   int maxiter;           // ������������ ���������� ��������
   double eps;            // �������� ��������� ������������� �������

   vector<double> b;      // ������ ������ �����
   vector<double> t;      // ��������������� ������ ��� ���
   vector<double> tt;     // ��������������� ������ ��� ���
   vector<double> rk1;    // ������ ������� �� ����. �������� ���
   vector<double> zk1;    // ������ ������ �� ����. �������� ���
   vector<double> AtAzk1; // ��������������� ������ ��� ���

   SLAE(int size, int _maxiter, double _eps)
   {
      maxiter = _maxiter;
      eps = _eps;

      b.resize(size);

      t.resize(size);
      tt.resize(size);
      rk1.resize(size);
      zk1.resize(size);
      AtAzk1.resize(size);
   }

   SLAE()
   {

   }

   // ����� ����������� ����������, ���������� ���������� ��������
   int ConjGradMethod(vector<double>& xk1, vector<double>& res, Matrix& mat)
   {
      for(int i = 0; i < mat.size; i++)
         res[i] = xk1[i] = 0;

      mat.MatVecMult(xk1, t, mat.bot_tr, mat.top_tr);       // t = A * x0
      mat.MatVecMult(b - t, rk1, mat.top_tr, mat.bot_tr);  // r0 = AT(f - A * x0)
      zk1 = rk1;

      int k = 1;
      while(k < maxiter)
      {
         mat.MatVecMult(zk1, t, mat.bot_tr, mat.top_tr);    // t = A * zk-1

         mat.MatVecMult(t, AtAzk1, mat.top_tr, mat.bot_tr); // AtAzk1 = At * A * zk-1
         double ak = (rk1 * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         double bk = rk1 * rk1;
         rk1 = rk1 - ak * AtAzk1;
         bk = (rk1 * rk1) / bk;
         zk1 = rk1 + bk * zk1;

         double disc = Norm(rk1) / Norm(b); // ������������� �������

         if(disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }
};
