Basis sets optimized for H2O molecule

Method: contraint ABS-0 by default, no orbital type constraint imposed

output is x_abs_y_0 obtained by reducing the x basis to size y 

Results generated with
======================

6-31g (ref size 9), reduced to 8, 7, 5:
python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 8
python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 7
python3 main.py --coord h2o --out out/test_3 --AO 6-31g --M_target 5

cc-pvdz (ref size 12), reduced to 10, *9*, 7, 5: 
python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 10
python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 9
python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 7
python3 main.py --coord h2o --out out/test_3 --AO cc-pvdz --M_target 5

cc-pvtz (ref size 22), reduced to 20, 15, *12*, 10:
python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 20
python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 10
python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 12
python3 main.py --coord h2o --out out/test_3 --AO cc-pvtz --M_target 15

cc-pvqz (ref size 35): Killed (basis too large)
python3 main.py --coord h2o --out out/test_3 --AO cc-pvqz --M_target 15

