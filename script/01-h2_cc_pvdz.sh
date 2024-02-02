#!/bin/bash

COORD=h2
AO=cc-pvdz

M_TARGETS="3 5 7"
OPTIONS="0 1 2"

for OPTION in $OPTIONS;
do
    for M_TARGET in $M_TARGETS;
    do
        python3 main.py --coord $COORD --AO $AO --M $M_TARGET --option $OPTION
    done
done
