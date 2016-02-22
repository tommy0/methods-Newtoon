#ifndef FUNCTION_H
#define FUNCTION_H
#include <vector>

class Function
{
    public:
        Function(double (*f)(std::vector<double>));
        double getValue(std::vector<double> inX);
    private:
        double (*f)(std::vector<double>);
};

#endif // FUNCTION_H
