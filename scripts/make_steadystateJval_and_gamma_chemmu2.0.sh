#!/bin/bash

savedirname='ce_elife_revision'
mkdir $savedirname

./plot_steadystateJval_and_gamma_v2.py ${savedirname}/ce_chemmu2.pdf 10 "$(printf %s, {10,50,100}000 )" ../../results/ce_elife_revision/ce_nof_noevr_s{10,50,100}k*chemmu2.0/da*


