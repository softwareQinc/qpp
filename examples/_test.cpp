// Used for testing, do not use it as an example
#include <iostream>
#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

int main()
{
    std::cout << "Testing...\n";
    PRINTLN("Testing debug messages");
    ERRORLN("Testing debug error messages");

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
}
