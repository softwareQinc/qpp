#!/bin/sh

find include examples unit_tests/tests -name '*.cpp' -o -name '*.c' -o -name '*.h'\
     -o -name '*.hpp' | xargs clang-format-3.9 -style=file -i
