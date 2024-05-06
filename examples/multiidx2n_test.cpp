#include <iostream>
#include "qpp/functions.hpp"
#include "qpp/internal/util.hpp"

int main() {

    constexpr int ndims = 3;
    int midx[ndims] = {1, 0, 1};
    int dims[ndims] = {2, 2, 2};

    std::cout << qpp::internal::multiidx2n(midx, ndims, dims) << std::endl;
}
