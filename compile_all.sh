#!/bin/bash

methods=$(grep -e '^\s*subroutine' harmonic.f90 | sed -e 's/(.*)//g' -e 's/  subroutine //g')
echo $methods
for method in $methods; do
  for n in 1 2 4; do
    for nrespa in `seq 55 5 100`; do
      make METHOD=$method NN=$n NRESPA=$nrespa clean
      make METHOD=$method NN=$n NRESPA=$nrespa
    done
  done
  mkdir -p $method
  mv harmonic-${method}* ${method}/
done
