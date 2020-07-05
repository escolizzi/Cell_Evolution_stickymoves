#!/bin/bash

# g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu0.1_*/data* 
# g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu0.5*/da* 
# g ../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_50c_f500_{1..3}.txt 
# g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu2.0*/da* 
# g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu5.0_*/da*

mkdir bc_1c_vs_50c_chemmu

chemmu='0.1'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf 250 0 g ../../results/bc_elife_revision/bc_g6_1c_f500_chemmu${chemmu}*/data* g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu${chemmu}*/da*
echo Saving file: bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf

chemmu='0.5'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf 250 0 g ../../results/bc_elife_revision/bc_g6_1c_f500_chemmu${chemmu}*/data* g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu${chemmu}*/da*
echo Saving file: bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf

chemmu='1.0'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf 250 0 g ../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_1c_f500_{1..10}.txt g ../../results/bc_g6_1c_50c_f500_onlydata/data_bc_g6_50c_f500_{1..10}.txt
echo Saving file: bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf

chemmu='2.0'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf 250 0 g ../../results/bc_elife_revision/bc_g6_1c_f500_chemmu${chemmu}*/data* g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu${chemmu}*/da*
echo Saving file: bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf

chemmu='5.0'
./plot_displacement_in_time_v4.py bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf 250 0 g ../../results/bc_elife_revision/bc_g6_1c_f500_chemmu${chemmu}*/data* g ../../results/bc_elife_revision/bc_g6_50c_f500_chemmu${chemmu}*/da*
echo Saving file: bc_1c_vs_50c_chemmu/bc_1c_vs_50c_chemmu${chemmu}.pdf


