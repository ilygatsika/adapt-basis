Basis sets optimized for H2O molecule

Method: contraint ABS-0 by default, no orbital type constraint imposed

output is x_abs_y_0 obtained by reducing the x basis to size y 

Results generated with script script/03-h2o.sh. In more detail:
==============================================================

6-31g (ref size 9), reduced to 8, 7, 5:

    python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 8
    python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 7
    python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 5

cc-pvdz (ref size 12), reduced to 10, *9*, 7, 5: 

    python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 10
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 9
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 7
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 5
    
cc-pvtz (ref size 22), reduced to 20, 15, *12*, 10, *9*:

    python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 20
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 15
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 12
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 10
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 9
    
Note: simulations below need cluster due to large memory allocation

cc-pvqz (ref size 35), reduced to 32, 30, 25, *22, 12, 9*
Note: 10-20 GB memory allocated

    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 32
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 30
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 25
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 22
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 12
    python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 9

cc-pv5z (ref size 51), reduced to 48, 40, 38, *35, 22, 12, 9*
Note: 40-55 GM memory allocated

    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 48
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 40
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 38
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 35
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 22
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 12
    python3 main.py --coord h2o --out out/test_3 --AO cc-pv5z --M_target 9

