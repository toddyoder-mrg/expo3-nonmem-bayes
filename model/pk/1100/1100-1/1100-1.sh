#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1100/1100-1

/opt/NONMEM/nm75/run/nmfe75 1100-1.ctl  1100-1.lst  -parafile=1100-1.pnm -maxlim=2
