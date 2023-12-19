// Statistics
// Source: ./examples/statistics.cpp

#include <iostream>
#include <vector>

#include "qpp/qpp.h"

int main() {
    using namespace qpp;

    // random variables
    std::vector<int> X{1, 2};
    std::vector<int> Y{1, 2, 3};

    // joint probability distribution
    rmat probXY(2, 3);
    probXY << 0.25, 0.25, 0, 0, 0.25, 0.25;
    std::vector<realT> probX = marginalX(probXY);
    std::vector<realT> probY = marginalY(probXY);

    std::cout << ">> ProbX: "
              << disp(marginalX(probXY), IOManipContainerOpts{}.set_sep(", "))
              << '\n';
    std::cout << ">> ProbY: "
              << disp(marginalY(probXY), IOManipContainerOpts{}.set_sep(", "))
              << '\n';

    std::cout << "Mean (X/Y): " << avg(probX, X) << " " << avg(probY, Y)
              << '\n';
    std::cout << "Standard deviation (X/Y): " << sigma(probX, X) << " "
              << sigma(probY, Y) << '\n';
    std::cout << "Variance (X/Y): " << var(probX, X) << " " << var(probY, Y)
              << '\n';
    std::cout << "Covariance: " << cov(probXY, X, Y) << '\n';

    // display a uniform probability distribution
    std::cout << "Uniform(5): "
              << disp(uniform(5), IOManipContainerOpts{}.set_sep(", ")) << '\n';
}
