#!/bin/bash

methods=$(cat methods.txt)
methods=xi_respa_NHC
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
