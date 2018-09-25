# PS regime not really tested otherwise
#!/bin/bash

CODEDIR=~/NEOART_2.0/

echo $CODEDIR

cd $CODEDIR

cp tests/input/input_test2.dat input.dat
echo 'COLLISIONAL PS REGIME:'

./neoart.x
echo 'ION HEAT FLUX SHOULD BE SHOULD BE -0.80'
