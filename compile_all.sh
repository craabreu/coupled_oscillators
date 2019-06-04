#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 NN NRESPA"
    exit 0
fi

methods=$(grep -e '^\s*subroutine' harmonic.f90 | sed -e 's/(.*)//g' -e 's/  subroutine //g')
echo $methods
for method in $methods; do
  make clean
  make METHOD=$method NN=$1 NRESPA=$2
done
