#!/bin/bash

# CHANGED TO FULL SCHEDULE TO AVOID TOO MANY THING RUNNING AT THE SAME TIME

FILE=../results/test_long_term_rk.txt
if test -f "$FILE"; then
    rm "$FILE"
fi

# THEN LONG TERM
N_S=6
NUMBER_OF_DAYS=100
TRUE_SYSTEM=4

########################
MODEL_SYSTEM=1 # SYSTEM

echo "Model 1" >> "$FILE"
echo "---------------------" >> "$FILE"
echo "Slope: 0.2" >> "$FILE"
price_slope=0.2
regime_slope=0.2
python simulate_mpc_rk_reuse.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope

# echo "Slope: 0.5" >> "$FILE"
# price_slope=0.5
# regime_slope=0.5
# python simulate_mpc_rk_reuse.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope

########################
# MODEL_SYSTEM=2 # SYSTEM
# echo " " >> "$FILE"
# echo "Model 2" >> "$FILE"
# echo "---------------------" >> "$FILE"
# echo "Slope: 0.2" >> "$FILE"
# price_slope=0.2
# regime_slope=0.2
# python simulate_mpc_rk.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope

# echo "Slope: 0.5" >> "$FILE"
# price_slope=0.5
# regime_slope=0.5
# python simulate_mpc_rk.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope

########################
# echo " " >> "$FILE"
# echo "Model 3" >> "$FILE"
# echo "---------------------" >> "$FILE"
# echo "Slope: 0.2" >> "$FILE"
# MODEL_SYSTEM=3 # SYSTEM
# price_slope=0.2
# regime_slope=0.2
# python simulate_mpc_rk_reuse.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope

# echo "Slope: 0.5" >> "$FILE"
# price_slope=0.5
# regime_slope=0.5
# python simulate_mpc_rk_reuse.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope


# LEFT TERMINAL