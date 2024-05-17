#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1100/1100-4

/opt/NONMEM/nm75/run/nmfe75 1100-4.ctl  1100-4.lst  -parafile=1100-4.pnm -maxlim=2
