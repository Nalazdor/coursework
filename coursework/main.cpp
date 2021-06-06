#include <iostream>
#include "HyperbolicalProblem.h"

int main()
{
   HyperbolicalProblem hp;

   hp.ReadFormGrid("grid.txt");
   hp.ReadFormTimeGrid("time_grid.txt");
   hp.InitializeMemory();
   hp.test = Test(6);

   hp.FormPortrait(hp.A);
   hp.FormPortrait(hp.M);
   hp.FormPortrait(hp.G);

   ofstream fout("result.txt");
   hp.IterateTime(fout);
   fout.close();
}