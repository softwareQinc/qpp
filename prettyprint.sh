#!/bin/sh

# $@ - List of directories

# Code beautifier with clang-format
# Recursively parses the directories passed as command line arguments

if test -z "$CLANG_FORMAT"; then
    echo "Please set the CLANG_FORMAT environment variable to point to the \
location of clang-format"
    exit 1
else
    if ! [ -x "$(command -v "$CLANG_FORMAT")" ]; then
        echo "Error: $CLANG_FORMAT executable not found." >&2
        exit 1
    fi
    echo "Code formatting with '$CLANG_FORMAT' the directories:"
fi

for directory in "$@"; do
    echo "$directory"
    find "$directory" \( -iname '*.cpp' -o -iname '*.c' -o -iname '*.h' \
        -o -iname '*.hpp' \) -exec "$CLANG_FORMAT" -style=file -i {} +
done
