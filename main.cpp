#include <iostream>
#include "Function.h"
#include "System.h"
#include "Input.h"

///TODO multi var function and turn X

int main()
{
    Function first(f1);
    Function second(f2);
    std::vector<Function> functions ={first, second};
    std::vector<double> x ={1,2};
    double epsilon=0.001;
    System equation(functions,x,epsilon);
    equation.solveSystem();
    x=equation.getX();
    for(unsigned int i=0;i<x.size();++i)
    {
        std::cout<<x[i]<<" ";
    }
    return 0;
}

