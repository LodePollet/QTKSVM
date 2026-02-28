#include <iostream>
#include <string>
#include <random>
#include <array>
#include <cmath>

int main()
{
    constexpr int dim = 3;
    int L = 9;
    int lin_l = 3;
    int n_basis = 3;
    std::array<size_t, dim + 1> plengths;
    
    plengths.front() = 1;
    std::fill(plengths.begin() + 1, plengths.end(), L);
    std::partial_sum(plengths.begin(), plengths.end(), plengths.begin(),
                        std::multiplies<>{});

    auto index = [&](size_t i) {
        size_t const pl_mcell = n_basis * pow(lin_l, 2.);
        size_t const rw_mcell = n_basis * lin_l;
        size_t const pl_lat = n_basis * pow(plengths[1], 2.);  
        size_t const rw_lat = n_basis * plengths[1];
        return (i/pl_mcell)*pl_lat + ((i%pl_mcell)/rw_mcell)*rw_lat +
                       (i % pl_mcell) % rw_mcell;
    };

    for (size_t i = 0; i < n_basis * pow(lin_l, dim); i++) {
        std::cout << i << " -> " << index(i) << std::endl;
    }
    return 0;
}
