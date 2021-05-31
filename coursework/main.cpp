#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   HyperbolicalProblem hp;
   hp.ReadFormGrid("grid.txt");
   hp.ReadFormTimeGrid("time_grid.txt");
}