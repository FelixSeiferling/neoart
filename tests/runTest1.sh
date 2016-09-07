#!/bin/bash

CODEDIR=~/NEOART_2.0/

cd $CODEDIR



echo 'ION TEMPERATURE GRADIENT:'
sed -e "s/dt = .0         ! temperature gradient/dt = 1.0         ! temperature gradient/" tests/input/input_test1.dat > input.dat
./neoart.x
echo 'HEAT FLUX SHOULD BE -1.80 / -0.68'
echo 'BOOTSTRAP SHOULD BE 2.86'
echo 'ION ROTATION SHOULD BE 1.17'
echo ' '

# echo 'ELECTRON TEMPERATURE GRADIENT:'
# sed -e "s/dt = .0         ! electron temperature gradient/dt = 1.0         ! electron temperature gradient/" tests/input/input_test1.dat > input.dat
# ./neoart.x
# echo 'HEAT FLUX SHOULD BE -3.34 / 0.0'
# echo 'BOOTSTRAP SHOULD BE 1.75'
# echo 'ION ROTATION SHOULD BE 0'
# echo ' '
# 
# echo 'ION PRESSURE GRADIENT:'
# sed -e "s/dp = .0         ! pressure gradient/dp = 1.0         ! pressure gradient/" tests/input/input_test1.dat > input.dat
# ./neoart.x
# echo 'HEAT FLUX SHOULD BE 1.53 / ~ 0'
# echo 'BOOTSTRAP SHOULD BE -2.44'
# echo 'ION ROTATION SHOULD BE 0'
# echo ' '
# 
# echo 'ELECTRON PRESSURE GRADIENT:'
# sed -e "s/dp = .0         ! electron pressure gradient/dp = 1.0         ! electron pressure gradient/" tests/input/input_test1.dat > input.dat
# ./neoart.x
# echo 'HEAT FLUX SHOULD BE 1.53 / 0.0'
# echo 'BOOTSTRAP SHOULD BE -2.44'
# echo 'ION ROTATION SHOULD BE 0'
# echo ' '

