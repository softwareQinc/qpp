// Used for testing, do not use it as an example
#include <iostream>
#include <fstream>
#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

int main()
{
    std::cout << "Testing...\n";

    std::cout << "Playing with exceptions...\n";
    try
    {
        throw exception::UndefinedType{__FILE__ + std::string(" at line ") +
                                       std::to_string(__LINE__)};
    }
    catch (exception::UndefinedType& e)
    {
        std::cout << e.what() << '\n';
        std::cout << "Exception type: " << e.type_description() << '\n';
    };


    for (idx i = 0; i < 10; ++i)
    {
        std::cout << rand(bigint(0), 10) << " ";
    }
    std::cout << '\n';

    std::ifstream fout("prng_state.txt");
    rdevs.load(fout);
    for (idx i = 0; i < 10; ++i)
    {
        std::cout << rand((bigint) 0, 100) << " ";
    }
}
