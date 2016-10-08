// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>

using namespace qpp;

int main()
{
    std::cout << "Testing..." << std::endl;
#define DEBUG
    ERRORLN("Oops");
#undef DEBUG
    ERRORLN("Not printed");
}
