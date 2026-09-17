#!/bin/bash

find . -name "*.c" | xargs clang-format -i
find . -name "*.h" | xargs clang-format -i
find . -name "*.cpp" | xargs clang-format -i
