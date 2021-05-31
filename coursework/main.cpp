#include <iostream>
#include "HyperbolicalProblem.h"

int main()
{
   HyperbolicalProblem hp;

   hp.ReadFormGrid("grid.txt");
   hp.ReadFormTimeGrid("time_grid.txt");
   hp.InitializeMemory();
   hp.test = Test(4);

   hp.FormPortrait(hp.A);
   hp.FormPortrait(hp.M);
   hp.FormPortrait(hp.G);

   hp.FillMatrices();
   hp.FillB(0);

   hp.AssembleGlobalMatrix();
   hp.AccountBound(0);

   hp.Solve();

   ofstream fout("result.txt");
   hp.PrintSolution(fout, 0);
   fout.close();
}