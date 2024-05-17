#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000/1000-1

/opt/NONMEM/nm75/run/nmfe75 1000-1.ctl  1000-1.lst  -parafile=1000-1.pnm -maxlim=2
