#!/bin/sh
# Code beautifier with clang-format
# Recursively parses the folder(s) passed as command line argument(s)

CLANG_FORMAT=clang-format-3.9

echo "Code formatting..."
for folder in "$@"
do
    echo $folder 
    find $folder \( -iname '*.cpp' -o -iname '*.c' -o -iname '*.h'\
         -o -iname '*.hpp' \) | xargs $CLANG_FORMAT -style=file -i
done
