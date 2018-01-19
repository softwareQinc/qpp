#!/bin/sh
# Code beautifier with clang-format

CLANG_FORMAT=clang-format-3.9

for folders in "$@"
do
    find $folders \( -iname '*.cpp' -o -iname '*.c' -o -iname '*.h'\
         -o -iname '*.hpp' \) | xargs $CLANG_FORMAT -style=file -i
done
