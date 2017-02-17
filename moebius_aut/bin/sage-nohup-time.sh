#! /bin/bash

#
# This script can be used to run this package on a remote server,
# via an ssh session. the standard/error output is written
# in the files err/out respectively.
#

export PYTHONPATH=$PYTHONPATH:../src/
export OUTPUT_PATH=/home/LOCAL/nlubbes/OUTPUT/

rm err out
nohup time nice /home/software/sage-7.3/sage -python ../src/moebius_aut/__main__.py  > out 2> err < /dev/null &
cat err
cat out

