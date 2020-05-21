#!/bin/bash


#N_S=6
NUMBER_OF_DAYS=20
TRUE_SYSTEM=3
price_slope=0.2
regime_slope=0.2

MODEL_SYSTEM=2
for N_S in 4 6 8 10 12 14 16 18 20
do
    python simulate_mpc_pure_dap.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope
done

MODEL_SYSTEM=3
for N_S in 4 6 8 10 12 14 16 18 20
do
    python simulate_mpc_pure_dap.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope
done