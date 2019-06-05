#!/bin/bash

methods=$(grep -e '^\s*subroutine' harmonic.f90 | sed -e 's/(.*)//g' -e 's/  subroutine //g')
echo $methods
for n in 1 2 4; do
  for nrespa in `seq 5 5 50`
    for method in $methods; do
      make METHOD=$method NN=$n NRESPA=$nrespa clean
      make METHOD=$method NN=$n NRESPA=$nrespa
      mkdir -p $method
      mv harmonic-${method}* ${method}/
    done
done
