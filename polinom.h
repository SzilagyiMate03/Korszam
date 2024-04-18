#include <vector>

template<typename T>
struct polinom
{
    std::vector<T> p;
    polinom<T>(std::vector<T> p_init) : p{p_init}{};
    polinom<T>() : p{}{};
    polinom<T>(size_t n) : p{std::vector<T> (n)}{};
};

template<typename T>
std::ostream& operator <<( std::ostream& o, polinom<T> const& a);
template<typename T>
polinom<T> operator+(polinom<T> const& a, polinom<T> const& b);
template<typename T>
polinom<T> operator-(polinom<T> const& a, polinom<T> const& b);
template<typename T>
polinom<T> dot(polinom<T> const& a, polinom<T> const& b);
template<typename T>
polinom<T> derivative(polinom<T> const& a, int order=1);
template<typename T>
T kiertekel(polinom<T> const& a, T x);
template<typename T>
T integrate(polinom<T> const& a, T begin, T end);

