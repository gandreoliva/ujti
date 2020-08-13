#!/bin/bash
doc="
--- -- -- -- ---._      UJTI - general geodesic ray-tracer
                o \`.                                   o     _
                    \          Cinespa UCR 2019  *  ( ) GANDREOLIVA.org
=======================================================================
Usage: ./make.sh task

tasks:
  -symb metric_name : generates .o, .mod objects for metric_name. The source files have to be located at gensrc/
  -num solver_src : generates .o, .mod objects for the solver in the specified source file solver_src

output:
  all outputs are located in bin/so/
8< - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - >8
"

case $1 in
-symb)
  if [ -z $2 ]; then echo "$doc"; echo "**Error: No metric specified"; exit 1; fi
  gfortran -c symbolic/geq-sphax-base.f90 -o bin/so/geq-sphax-base.o -Jbin/so/
  gfortran -c gensrc/$2.f90 -o bin/so/$2.o -Jbin/so/
  ;;
-num)
  if [ -z $2 ]; then echo "$doc"; echo "**Error: no solver source file specified"; exit 1; fi
  fnameext=$(basename $2)
  gfortran -c $2 -o bin/so/"${fnameext%.*}.o" -Jbin/so/
  ;;

*) echo "$doc"
  ;;
esac
