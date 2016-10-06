// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>

using namespace qpp;

int main()
{
    std::cout << "Testing..." << std::endl;

    std::cout << std::boolalpha << isprime(10000000019) << std::endl;
}
