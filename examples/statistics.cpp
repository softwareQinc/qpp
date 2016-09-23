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
    dmat probXY(2, 3);
    probXY << 0.25, 0.25, 0, 0, 0.25, 0.25;
    std::vector<double> probX = marginalX(probXY);
    std::vector<double> probY = marginalY(probXY);

    cout << ">> ProbX: " << disp(marginalX(probXY), ", ") << endl;
    cout << ">> ProbY: " << disp(marginalY(probXY), ", ") << endl;

    cout << "Mean (X/Y): " << avg(probX, X) << " " << avg(probY, Y) << endl;
    cout << "Standard deviation (X/Y):" << sigma(probX, X) << " "
         << sigma(probY, Y) << endl;
    cout << "Variance (X/Y): " << var(probX, X) << " " << var(probY, Y) << endl;
    cout << "Covariance: " << cov(probXY, X, Y) << endl;

    // display a uniform probability distribution
    cout << "Uniform(5): " << disp(uniform(5), ", ") << endl;
}
