#!/bin/bash
gprof ../build/FEM_PROFILE ../build/gmon.out > ../build/analysis.txt
gprof2dot ../build/analysis.txt -s > ../build/profiling.dot
xdot ../build/profiling.dot
rm -f ../build/gmon.out
rm -f ../build/analysis.txt
