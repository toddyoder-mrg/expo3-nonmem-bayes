#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1100/1100-3

/opt/NONMEM/nm75/run/nmfe75 1100-3.ctl  1100-3.lst  -parafile=1100-3.pnm -maxlim=2
