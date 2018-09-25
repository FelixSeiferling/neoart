#!/bin/bash
#this should work modify directories accordingly
CODEDIR=~/NEOART_2.0/

cd $CODEDIR



echo 'ION TEMPERATURE GRADIENT:'
sed -e "s/dt = .0         ! temperature gradient/dt = 1.0         ! temperature gradient/" tests/input/input_test1.dat > input.dat
./neoart.x
echo 'HEAT FLUX SHOULD BE -1.819 / -0.6816'
echo 'BOOTSTRAP SHOULD BE 2.777'
echo 'ION ROTATION SHOULD BE 1.164'
echo ' '

echo 'ELECTRON TEMPERATURE GRADIENT:'
sed -e "s/dt = .0         ! electron temperature gradient/dt = 1.0         ! electron temperature gradient/" tests/input/input_test1.dat > input.dat
./neoart.x
echo 'HEAT FLUX SHOULD BE -3.422 / 0.0'
echo 'BOOTSTRAP SHOULD BE 1.709'
echo 'ION ROTATION SHOULD BE 0'
echo ' '

echo 'ION PRESSURE GRADIENT:'
sed -e "s/dp = .0         ! pressure gradient/dp = 1.0         ! pressure gradient/" tests/input/input_test1.dat > input.dat
./neoart.x
echo 'HEAT FLUX SHOULD BE 1.562 / 0'
echo 'BOOTSTRAP SHOULD BE -2.384'
echo 'ION ROTATION SHOULD BE 0'
echo ' '

echo 'ELECTRON PRESSURE GRADIENT:'
sed -e "s/dp = .0         ! electron pressure gradient/dp = 1.0         ! electron pressure gradient/" tests/input/input_test1.dat > input.dat
./neoart.x
echo 'HEAT FLUX SHOULD BE 1.562 / 0.0'
echo 'BOOTSTRAP SHOULD BE -2.384'
echo 'ION ROTATION SHOULD BE 0'
echo ' '

