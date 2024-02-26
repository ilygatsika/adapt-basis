#!/bin/bash

# run in background with 
#   >>> nohup ./script/03-h2o.sh > out/test_3/out.log &
#
# Execution on server recommended (150 GB allocated memory)

COORD=h2o
OUT_DIR=out/test_3

DMS="HF CISD"
IS_ATOM="0 1"
M_TARGETS="7 8 9 10 11 12 13"
AOS='cc-pvdz cc-pvtz cc-pvqz cc-pv5z cc-pv6z'

for AO in $AOS;
do
    for M_TARGET in $M_TARGETS;
    do
        for OPTION in $IS_ATOM;
        do
            for DM in $DMS;
            do
                python3 main.py --coord $COORD --dm $DM --out $OUT_DIR --AO $AO --M_target $M_TARGET --gram_atom $OPTION
            done
        done
    done
done

