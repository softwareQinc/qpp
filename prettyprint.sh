#!/bin/sh

# Code beautifier with clang-format
# Recursively parses the directories passed as command line arguments

# Arguments:
#
# $@ - List of directories

# Exit immediately if a command exits with a non-zero status
set -e

# --- Configuration ---
# Allow user to override the command via env var, otherwise default to 'clang-format'
FORMATTER="${CLANG_FORMAT:-clang-format}"

# --- Validation ---
# Check if the formatter exists in the system path.
if ! command -v "$FORMATTER" >/dev/null 2>&1; then
    echo "Error: '$FORMATTER' not found." >&2
    echo "Please install clang-format or set the CLANG_FORMAT environment variable." >&2
    exit 1
fi

# Check if arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 [directory1] [directory2] ..."
    exit 1
fi

# --- Execution ---
echo "Formatting code using: $(command -v "$FORMATTER")"

for directory in "$@"; do
    if [ -d "$directory" ]; then
        echo "Processing: $directory"

        # We use a single 'find' command to handle everything.
        # 1. -name .git -prune : Ignores the .git directory entirely (efficiency).
        # 2. -type f : Only looks for files.
        # 3. \( ... \) : Groups the file extensions.
        # 4. -exec ... + : Passes multiple files to clang-format at once (batching).

        find "$directory" \
            -name ".git" -prune -o \
            -type f \( \
            -name '*.cpp' -o \
            -name '*.c++' -o \
            -name '*.cxx' -o \
            -name '*.cc' -o \
            -name '*.c' -o \
            -name '*.h' -o \
            -name '*.hpp' -o \
            -name '*.hh' -o \
            -name '*.hxx' \
            \) -exec "$FORMATTER" -style=file -i {} +
    else
        echo "Warning: '$directory' is not a valid directory. Skipping." >&2
    fi
done

echo "Done."
