#!/bin/bash

Np=(2 4 8 16 32 64 128 256 512 1024 2048)
Ne=(1000000 1000000 1000000 1000000 100000 100000 100000 10000 10000 1000 1000)

for i in {0..10}; do
  printf "${Np[$i]}: "
  ./bin/benchmark_fj antikt 0.4 ${Np[$i]} ${Ne[$i]} | tail -1
done
