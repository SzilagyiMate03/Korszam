#include <vector>

template<typename T>
struct polinom
{
    std::vector<T> p;
    polinom(std::vector<T> p_init) : p{p_init}{};
};