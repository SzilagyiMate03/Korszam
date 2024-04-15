#include <iostream>
#include <vector>
#include "polinom.h"

int main()
{
    polinom<int> a = std::vector<int> {1,2,3};
    std::cout << a.p[0] << std::endl;
    return 0;
}