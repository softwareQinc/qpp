#!/bin/sh
# Code beautifier with clang-format
# Recursively parses the folder(s) passed as command line argument(s)

if test -z "$CLANG_FORMAT"; then
    echo "Please set the CLANG_FORMAT environment variable to point to the \
location of clang-format"
    exit 1
else
    if ! [ -x "$(command -v "$CLANG_FORMAT")" ]; then
        echo "Error: $CLANG_FORMAT executable not found." >&2
        exit 1
    fi
    echo "Code formatting with '$CLANG_FORMAT' the folders:"
fi

for folder in "$@"; do
    echo "$folder"
    find "$folder" \( -iname '*.cpp' -o -iname '*.c' -o -iname '*.h'\
         -o -iname '*.hpp' \) -exec "$CLANG_FORMAT" -style=file -i {} +
done
