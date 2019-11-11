#!/bin/sh

cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -Wno-dev -S . -B build
