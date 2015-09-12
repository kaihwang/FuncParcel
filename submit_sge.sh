#!/bin/bash

# bash script to use sge to submit jobs iterating modularity detection

echo 'python /home/despoB/kaihwang/bin/FuncParcel/create_avemat.py' > ~/tmp/avg_partition.sh
qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/avg_partition.sh

# for Subject in 1103 1220 1306 1223 1314 1311 1318 1313 1326 1325 1328 1329 1333 1331 1335 1338 1336 1339 1337 1344 1340 128 162 163 168 176 b116 b117 b120 b121 b122 b138 b143 b153; do

# 	sed "s/128/${Subject}/g" < python_brainx_mod.sh > ~/tmp/python_brainx_${Subject}.sh
# 	qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/python_brainx_${Subject}.sh

# done