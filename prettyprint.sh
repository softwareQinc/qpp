#!/bin/sh

find include examples -name '*.cpp' -o -name '*.c' -o -name '*.h' -o -name '*.hpp' \
| xargs clang-format-3.9 -style=file -i
