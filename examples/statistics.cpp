// Statistics
// Source: ./examples/statistics.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    // random variables 
    std::vector<int> X{1, 2};
    std::vector<int> Y{1, 2, 3};

    // joint probability distribution
    dmat probXY(2,3);
    probXY << 0.25, 0.25, 0, 0, 0.25, 0.25;
    std::vector<double> probX = marginalX(probXY);
    std::vector<double> probY = marginalY(probXY);
    
    std::cout << ">> ProbX: " << disp(marginalX(probXY),", ") << std::endl;
    std::cout << ">> ProbY: " << disp(marginalY(probXY),", ") << std::endl;

    std::cout << "Mean (X/Y): " << avg(probX, X) << " " 
            << avg(probY, Y) << std::endl;
    std::cout << "Standard deviation (X/Y):" << sigma(probX, X) << " " 
            << sigma(probY, Y) << std::endl;
    std::cout << "Variance (X/Y): " << var(probX, X) << " " 
            << var(probY, Y) << std::endl;
    std::cout << "Covariance: " << cov(probXY, X, Y) << std::endl;

    // display a uniform probability distribution
    std::cout << "Uniform(5): " << disp(uniform(5),", ") << std::endl;
}
