// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>

using namespace qpp;

int main()
{
    std::cout << "Testing..." << std::endl;
    std::cout << std::numeric_limits<bigint>::min() << std::endl;
    std::cout << std::numeric_limits<bigint>::max() << std::endl;
    std::cout << 10 % 3 << std::endl;
}
