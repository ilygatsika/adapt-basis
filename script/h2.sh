#!/bin/bash

# run in background with
#   >>> nohup ./script/h2.sh > out/h2/out.log &
 
COORD=h2
OUT_DIR=out/h2

AO='cc-pvdz'
DM="H2_alpha_rdm"

IS_ATOM="0 1"
M_TARGETS="2 3 4 5 6"

for M_TARGET in $M_TARGETS;
do
    for OPTION in $IS_ATOM;
    do
        python3 main.py --coord $COORD --dm $DM --out $OUT_DIR --AO $AO --M_target $M_TARGET --gram_atom $OPTION
    done
done

DMS="HF CISD"
IS_ATOM="0 1"
M_TARGETS="2 3 4 5 6"
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

