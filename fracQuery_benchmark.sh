#!/bin/bash

seed=$1
for method in CCA liger conos
  do
  for frac in $(seq 0.1 0.1 1)
    do
    Rscript ~/multiOmic_benchmark/run_integration_fracQuery.R ~/my_data/10X_data/PBMC_SCElist_20191105.RDS $method $frac $seed
    done
  done