#!/bin/bash

#../../results/bc_elife_revision/bc_g6_50c_f500_mup0.1_*/da*
#../../results/bc_elife_revision/bc_g6_50c_f500_mup1.0_*/da*
#../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_50c_f500_*.txt
#../../results/bc_elife_revision/bc_g6_50c_f500_mup5.0_*/da*
#../../results/bc_elife_revision/bc_g6_50c_f500_mup10.0_*/da*


mkdir bc_1c_vs_50c_mup

for mup in 0.1 1.0 2.0 3.0 4.0 5.0 10.0;
do

./plot_displacement_in_time_v4.py bc_1c_vs_50c_mup/bc_1c_vs_50c_mup${mup}.pdf 250 0 g ../../results/bc_elife_revision/bc_g6_1c_f500_mup${mup}*/data* g ../../results/bc_elife_revision/bc_g6_50c_f500_mup${mup}*/da*
echo Saving file: bc_1c_vs_50c_mup/bc_1c_vs_50c_mup${mup}.pdf

done 

mup='3.0'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_mup/bc_1c_vs_50c_default_mup${mup}.pdf 250 0 g ../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_1c_f500_{1..10}.txt g ../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_50c_f500_{1..10}.txt
echo Saving file: bc_1c_vs_50c_mup/bc_1c_vs_50c_default_mup${mup}.pdf
