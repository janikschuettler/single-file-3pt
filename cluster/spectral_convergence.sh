#!/bin/bash
#$ -S /bin/bash
#$ -q *
#$ -N conv_N=7
#$ -l h_rt=500:00:00
#$ -t 1-11:1
#$ -cwd

. /etc/profile.d/modules.sh
module add shared
module add gcc/9.3

# N=(6 7 8 9 10 6 7 8 9 10 6 7 8 9 10 6 7 8 9 10 6 7 8 9 10 )
# tag=(2)
# max_many_eig=(50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200)
max_many_eig=(100 110 120 130 140 150 160 170 180 190 200)

./conv 7 3 ${max_many_eig[$SGE_TASK_ID-1]} 0.5