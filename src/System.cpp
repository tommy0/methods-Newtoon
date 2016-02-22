#include "System.h"

System::System(std::vector<Function>arrF, std::vector<double> x, double epsilon)
{
    this->arrF=arrF;
    setX(x);
    setEpsilon(epsilon);
    setArrDerivative();
    setArrValue();
    epsilon=0;
    setDx();
}

void System::setArrDerivative()
{
    arrDerivative.resize(arrF.size());
    for(unsigned int i=0; i<arrF.size(); ++i)
    {
        arrDerivative[i].resize(arrF.size());
        for(unsigned int j=0; j<arrF.size(); ++j)
            arrDerivative[i][j]=getDerivative(i,j);
    }

}

double System::getDerivative(const unsigned int i,unsigned int j)
{
    std::vector<double> temp=x;
    for(unsigned k=0;k<temp.size();++k)
        temp[k]+=epsilon;
    double t=(arrF[i].getValue(temp)-arrF[i].getValue(x))/epsilon;
    int k;
      return t;
}

void System::setArrValue()
{
    arrValue.resize(arrF.size());
    for(unsigned int i=0; i<arrF.size(); ++i)
        arrValue[i]=arrF[i].getValue(x);
}

void System::setX(std::vector<double> x)
{
    this->x=x;
}

void System::setEpsilon(double epsilon)
{
    this->epsilon=epsilon;
}

void System::setDx()
{
    dx.resize(arrF.size());
    for(unsigned int i=0; i<arrF.size();++i)
    {
        dx[i]=0;
    }
}

double System::getDet(std::vector<std::vector<double>> matrix)
{
    int l;
    double d;
    double sum11=1,sum12=0, sum21=1, sum22=0;
    for (unsigned int i=0; i<matrix.size(); i++)
    {
        sum11=1;
        l=2*matrix.size()-1-i;
        sum21=1;
        for (unsigned int j=0; j<matrix.size(); j++)
        {
            sum21*=matrix[j][l%(matrix.size())];
            l--;
            sum11*=matrix[j][(j+i)%(matrix.size())];
        }
        sum22+=sum21;
        sum12+=sum11;
    }
    d=sum12-sum22;
    return d;
}

std::vector<double> System::getX()
{
    return x;
}

void System::solveKramerDx()
{
    double det=getDet(arrDerivative);
    for(unsigned int i=0; i<dx.size();++i)
    {
        std::vector<std::vector<double>> temp=arrDerivative;
        temp[i]=arrValue;
        dx[i]=getDet(temp)/det;
    }
}

bool System::isAbsDxLessEpsilon()
{
   std::vector<double>::const_iterator largest = max_element( dx.begin(), dx.end() );
   if(*largest<epsilon) return false;
   return true;
}

void System::solveSystem()
{
    do
    {
        for(unsigned int i=0; i<x.size();++i)
            x[i]+=dx[i];
        setArrDerivative();
        setArrValue();
        solveKramerDx();
    }
    while(isAbsDxLessEpsilon());
}
