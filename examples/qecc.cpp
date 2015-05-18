// Quantum error correcting codes
// Source: ./examples/qecc.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    ket a0 = codes.codeword(Codes::Type::FIVE_QUBIT, 0);
    ket a1 = codes.codeword(Codes::Type::FIVE_QUBIT, 1);

    ket b0 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 0);
    ket b1 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 1);

    ket c0 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 0);
    ket c1 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 1);

    cout << ">> [[5, 1, 3]] Five qubit code. ";
    cout << ">> Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(a0) * a1) << endl;

    cout << ">> [[7, 1, 3]] Seven qubit Steane code. ";
    cout << ">> Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(b0) * b1) << endl;

    cout << ">> [[9, 1, 3]] Nine qubit Shor code. ";
    cout << ">> Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(c0) * c1) << endl;
}
