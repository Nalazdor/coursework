#pragma once
#pragma once
using namespace std;

class Test
{
public:
   int test_n = 0;            // Номер теста
   double lambda = 1;         
   double sigma  = 0;
   double chi    = 0;

   Test(const int& t_N) : test_n(t_N) {};
   Test() { };

   double f(const double& x, const double& t)
   {
      return -1 * divlambdagrad(x, t) +
         sigma * dudt(x, t) + chi * d2udt2(x, t);
      //return -1 * divlambdagrad(x, t) + sigma * u(x, t);
   }

   // Точное решение
   double u(const double& x, const double& t)
   {
      switch(test_n)
      {
         case 0: return 2.0;
         case 1: return x;
         case 2: return x * x;
         case 3: return x * x * x;
         case 4: return x * x * x * x;

         case 5: return t;
         case 6: return t * t;
         case 7: return t * t * t;
         case 8: return t * t * t * t;
         
      };
   }

   double divlambdagrad(const double& x, const double& t)
   {

      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 2;
         case 3: return 6 * x;
         case 4: return 12 * x * x;

         case 5: return 0;
         case 6: return 0;
         case 7: return 0;
         case 8: return 0;
      };

   }

   // Производная точного решения по t
   double dudt(const double& x, const double& t)
   {
      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 0;
         case 3: return 0;
         case 4: return 0;
              
         case 5: return 1;
         case 6: return 2 * t;
         case 7: return 3 * t * t;
         case 8: return 4 * t * t * t;
      };
   }

   // Вторая производная точного решения по t
   double d2udt2(const double& x, const double& t)
   {
      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 0;
         case 3: return 0;
         case 4: return 0;
              
         case 5: return 0;
         case 6: return 2;
         case 7: return 6 * t;
         case 8: return 12 * t * t;
      };
   }
};