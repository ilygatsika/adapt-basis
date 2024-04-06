#!/bin/bash

# run in background with 
#
#   >>> nohup ./script/all.sh h2 > out/h2/out.log &
#   >>> nohup ./script/all.sh h2o > out/h2o/out.log &
#   >>> nohup ./script/all.sh fh > out/fh/out.log &
#   >>> nohup ./script/all.sh lih > out/lih/out.log &
#
# Execution on server recommended (150 GB allocated memory)

system=$1
COORD=$system
OUT_DIR=out/$system

DM="HF"
AOS='cc-pvdz cc-pvtz cc-pvqz cc-pv5z cc-pv6z'

for AO in $AOS;
do
    python3 main.py --coord $COORD --dm $DM --out $OUT_DIR --AO $AO
done

git add out/$system
git commit -m 'ljll server results $system'
git push


