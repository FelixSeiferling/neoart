#!/bin/bash

CODEDIR=~/NEOART_2.0/
order=3
cd $CODEDIR
cp ~/Simulations/S_32_new/vnlin_moments_AV.dat $CODEDIR/vnlin_moments_AV.dat

# echo 'DIIID PARAMTERS'
# echo 'VNLIN_OFF: '
# sed -e "s/vnlin_order = 0 !/vnlin_order = 0 !/" input_vnlin_DIIID_style.dat > input.dat
# ./neoart.x
# echo 'VNLIN_ON: '
# # for order in -1 -2 -3 ;do
# # echo 'ORDER ' $order 'only'
# # sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_ASDEX_style.dat > input.dat
# # ./neoart.x
# # done
# order=3
# echo 'up to order' $order 
# sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_DIIID_style.dat > input.dat
# ./neoart.x
# 
# 
# echo 'ASDEX PARAMTERS'
# echo 'VNLIN_OFF: '
# sed -e "s/vnlin_order = 0 !/vnlin_order = 0 !/" input_vnlin_ASDEX_style.dat > input.dat
# ./neoart.x
# echo 'VNLIN_ON: '
# # for order in -1 -2 -3 ;do
# # echo 'ORDER ' $order 'only'
# # sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_ASDEX_style.dat > input.dat
# # ./neoart.x
# # done
# order=3
# echo 'up to order' $order 
# sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_ASDEX_style.dat > input.dat
# ./neoart.x
# 

echo 'ITER PARAMTERS'
echo 'VNLIN_OFF: '
sed -e "s/vnlin_order = 0 !/vnlin_order = 0 !/" input_vnlin_ITER_style.dat > input.dat
./neoart.x
echo 'VNLIN_ON: '
# for order in -1 -2 -3 ;do
# echo 'ORDER ' $order 'only'
# sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_ITER_style.dat > input.dat
# ./neoart.x
# done
order=3
echo 'up to order' $order 
sed -e "s/vnlin_order = 0 !/vnlin_order = $order !/" input_vnlin_ITER_style.dat > input.dat
./neoart.x


