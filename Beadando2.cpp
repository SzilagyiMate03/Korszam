#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "polinom.h"

template<typename T>
std::ostream& operator <<( std::ostream& o, polinom<T> const& a)
{
    o << "[ ";
    for (T const& element : a.p) {o << element << " ";}
    o << "]";
    return o;
}


template<typename T>
polinom<T> operator+(polinom<T> const& a, polinom<T> const& b)
{
    size_t result_size = std::max(a.p.size(), b.p.size());
    polinom<T> result(result_size);

    for (size_t i = 0; i < result_size; ++i) {
        T a_coefficient = (i < a.p.size()) ? a.p[i] : 0;
        T b_coefficient = (i < b.p.size()) ? b.p[i] : 0;
        result.p[i] = a_coefficient + b_coefficient;
    }

    return result;
}

template<typename T>
polinom<T> operator-(polinom<T> const& a, polinom<T> const& b)
{
    size_t result_size = std::max(a.p.size(), b.p.size());
    polinom<T> result(result_size);

    for (size_t i = 0; i < result_size; ++i) {
        T a_coefficient = (i < a.p.size()) ? a.p[i] : 0;
        T b_coefficient = (i < b.p.size()) ? b.p[i] : 0;
        result.p[i] = a_coefficient - b_coefficient;
    }

    return result;
}

template<typename T>
polinom<T> dot(polinom<T> const& a, polinom<T> const& b)
{
    polinom<T> result(a.p.size()+b.p.size()-1);
    for (size_t i = 0; i < a.p.size(); i++)
    {
        for (size_t j = 0; j < b.p.size(); j++)
        {
            result.p[i+j] += a.p[i]*b.p[j];
        }
    }
    return result;
}

template<typename T>
polinom<T> derivative(polinom<T> const& a, int order=1)
{
    polinom<T> result(a.p.size()-order);
    for (int i = result.p.size()-1; i>=0; i--)
    {
        result.p[i] = a.p[i + order] * tgamma(i + order+1) / tgamma(i+1);
    }
    return result;
}

template<typename T>
T kiertekel(polinom<T> const& a, T x)
{
    T result = 0;
    for (int i=0; i<a.p.size();i++)
    {
        result+=pow(x,i)*a.p[i];
    }
    return result;
}

template<typename T>
T integrate(polinom<T> const& a, T begin, T end)
{
    T result = 0;
    for (size_t i =0; i <= a.p.size();i++)
    {
        result += a.p[i]*(pow(end,i+1)-pow(begin,i+1))/(i+1);
    }
    return result;
}

int main()
{
    polinom<int> a{std::vector<int> {1,2,3,4}};
    polinom<int> b{std::vector<int> {5,7,9}};
    polinom<int> null_polinom(4);
    polinom<int> c = b + null_polinom;
    polinom<int> d = a - b;
    polinom<int> e = dot(a,b);
    polinom<int> f = derivative(a,2);
    int g = kiertekel(a,2);
    int h = integrate(a,0,1);
    std::cout << c << std::endl;
    std::cout << d << std::endl;
    std::cout << e << std::endl;
    std::cout << f << std::endl;
    std::cout << g << std::endl;
    std::cout << h << std::endl;
    return 0;
}