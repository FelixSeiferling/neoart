#!/bin/bash

CODEDIR=~/NEOART_2.0/
order=3
cd $CODEDIR

echo 'RLT scan with ITER parameters:'
for RLT in 2 4 6 8 10; do
echo 'RLT='$RLT
echo 'VNLIN_OFF: '
dp=$(echo " $RLT+2.2" | bc)
sed -e "s/vnlin_order = 0 !/vnlin_order = 0 !/" -e "s/dt = 5.5         ! electron/dt = $RLT         ! electron !/" -e "s/dp = 7.7         ! electron/dp = $dp         ! electron/" -e "s/dt = 6.9         ! temperature/dt = $RLT         ! temperature/"  -e "s/dp = 9.1         ! pressure/dp = $dp         ! pressure/" input_vnlin_ITER_style.dat > input.dat
cp ~/Simulations/RLTscan/RLT_$RLT/vnlin_moments_AV.dat $CODEDIR/vnlin_moments_AV.dat
./neoart.x
echo 'VNLIN_ON: '
sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" -e "s/dt = 5.5         ! electron/dt = $RLT         ! electron !/" -e "s/dp = 7.7         ! electron/dp = $dp         ! electron/" -e "s/dt = 6.9         ! temperature/dt = $RLT         ! temperature/"  -e "s/dp = 9.1         ! pressure/dp = $dp         ! pressure/" input_vnlin_ITER_style.dat > input.dat
./neoart.x
done
