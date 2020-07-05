#!/bin/bash

savedirname='ce_elife_revision'
mkdir $savedirname

./plot_steadystateJval_and_gamma_v2.py ${savedirname}/ce_mup1.pdf 10 "$(printf %s, {10,50,100}000 )" ../../results/ce_elife_revision/ce_nof_noevr_s{10,50,100}k*mup1.0/da*

./plot_steadystateJval_and_gamma_v2.py ${savedirname}/ce_mup5.pdf 10 "$(printf %s, {10,50,100}000 )" ../../results/ce_elife_revision/ce_nof_noevr_s{10,50,100}k*mup5.0/da*

