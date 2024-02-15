#!/bin/bash

# run in background with ./script/03-h2o.sh &
# Attention should be run in a cluster (150 GB allocated memory)

COORD=h2o
OUT_DIR=out/test_3

AO=6-31g
M_TARGETS="8 7 5"
for M_TARGET in $M_TARGETS;
do
    python3 main.py --coord $COORD --out $OUT_DIR --AO $AO --M_target $M_TARGET
done

AO=cc-pvdz
M_TARGETS="10 9 7 5"
for M_TARGET in $M_TARGETS;
do
    python3 main.py --coord $COORD --out $OUT_DIR --AO $AO --M_target $M_TARGET
done

AO=cc-pvtz
M_TARGETS="20 15 12 10 9"
for M_TARGET in $M_TARGETS;
do
    python3 main.py --coord $COORD --out $OUT_DIR --AO $AO --M_target $M_TARGET
done

AO=cc-pvqz
M_TARGETS="32 30 25 22 12 9"
for M_TARGET in $M_TARGETS;
do
    python3 main.py --coord $COORD --out $OUT_DIR --AO $AO --M_target $M_TARGET
done

AO=cc-pv5z
M_TARGETS="48 40 38 35 22 12 9"
for M_TARGET in $M_TARGETS;
do
    python3 main.py --coord $COORD --out $OUT_DIR --AO $AO --M_target $M_TARGET
done

