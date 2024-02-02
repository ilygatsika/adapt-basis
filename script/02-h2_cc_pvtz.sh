#!/bin/bash

COORD=h2
AO=cc-pvtz

M_TARGETS="8 10 12"
OPTIONS="0 1 2"

for OPTION in $OPTIONS;
do
    for M_TARGET in $M_TARGETS;
    do
        python3 main.py --coord $COORD --AO $AO --M $M_TARGET --option $OPTION
    done
done
