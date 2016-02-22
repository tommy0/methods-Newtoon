#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <algorithm>
#include "Function.h"

class System
{
    public:
        System(std::vector<Function>arrF, std::vector<double> x, double epsilon);
        void getResult(std::vector<double> &res);
        std::vector<double> getX();
        void solveSystem();
        double getDerivative(const unsigned int i,const unsigned int j);
        double getDet(std::vector<std::vector<double>> matrix);
    private:
        double epsilon;
        std::vector<Function>arrF;
        std::vector<double> x;
        std::vector<std::vector<double>> arrDerivative;
        std::vector<double> dx;
        std::vector<double> arrValue;
        void setX(std::vector<double> x);
        void setEpsilon(double epsilon);
        void setArrDerivative();
        void setDx();
        void setArrValue();
        void solveKramerDx();
        bool isAbsDxLessEpsilon();
};

#endif // SYSTEM_H
