#!/bin/bash
FILES=$(find . -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.cc" -o -name "*.c" -o -name "*.h" \) -not -path "./lib/*" -not -path "./build/*")
echo "+ Running clang-format on the following files:"
echo $FILES
clang-format -i $FILES