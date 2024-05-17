#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1100/1100-2

/opt/NONMEM/nm75/run/nmfe75 1100-2.ctl  1100-2.lst  -parafile=1100-2.pnm -maxlim=2
