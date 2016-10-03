// Randomness
// Source: ./examples/randomness.cpp
#include <qpp.h>

using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    cout << ">> Generating a random ket on D = 5" << endl;
    ket rket = randket(5);
    cout << disp(rket) << endl;

    std::vector<double> probs = abssq(rket);
    cout << ">> Probabilities: " << disp(probs, ", ") << endl;

    cout << ">> Sum of the probabilities: ";
    cout << sum(probs.begin(), probs.end()) << endl;
}
