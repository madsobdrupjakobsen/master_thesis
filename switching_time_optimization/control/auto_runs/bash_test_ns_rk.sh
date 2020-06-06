#!/bin/bash

FILE=../results/test_ns_rk.txt
if test -f "$FILE"; then
    rm "$FILE"
fi

# NOT USED ANYMORE!!
N_S=6
NUMBER_OF_DAYS=20
TRUE_SYSTEM=3
price_slope=0.2
regime_slope=0.2

# echo "Model 1" >> "$FILE"
# echo "---------------------" >> "$FILE"
# MODEL_SYSTEM=1
# for N_S in 4 6 8 10 12 14 16 18 20
# do
#     echo "NS: $N_S" >> "$FILE"
#     python simulate_mpc_rk.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope
# done

echo " " >> "$FILE"
echo "Model 2" >> "$FILE"
echo "---------------------" >> "$FILE"
MODEL_SYSTEM=2
for N_S in 18
do
    echo "NS: $N_S" >> "$FILE"
    python simulate_mpc_rk.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope
done

# echo " " >> "$FILE"
# echo "Model 3" >> "$FILE"
# echo "---------------------" >> "$FILE"
# MODEL_SYSTEM=3
# for N_S in 4 6 8 10 12 14 16 18 20
# do
#     echo "NS: $N_S" >> "$FILE"
#     python simulate_mpc_rk.py $MODEL_SYSTEM $NUMBER_OF_DAYS $N_S $TRUE_SYSTEM $price_slope $regime_slope
# done


# LEFT terminal