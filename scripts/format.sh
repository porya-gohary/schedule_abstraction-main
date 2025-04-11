#!/bin/bash
find . -regex '.*\.\(cpp\|hpp\|cc\|c\|h\)' -not -path "./build/*" -exec clang-format -i {} +