#include "Function.h"

Function::Function(double (*f)(std::vector<double>))
{
    this->f=f;
}

double Function::getValue(std::vector<double> inX)
{
    return f(inX);
}
