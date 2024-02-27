#!/bin/bash

# run in background with
#   >>> nohup ./script/04-h2.sh > out/test_4/out.log &
 
COORD=fh
OUT_DIR=tmp/

DMS="HF"
IS_ATOM="1"
M_TARGETS="2 3 4 5"
AOS='cc-pvdz cc-pvtz cc-pvqz cc-pv5z'

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

