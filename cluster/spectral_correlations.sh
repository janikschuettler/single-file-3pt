#!/bin/bash
#$ -S /bin/bash
#$ -q *
#$ -N nm-t1,x0.0
#$ -l h_rt=80:00:00
#$ -t 1-21:1
#$ -cwd

. /etc/profile.d/modules.sh
module add shared
module add gcc/9.3

dir=simulation

n_threads=50
max_eig=400
nbins_t=100

N=5
tag=1
x=0.0
x1=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

./nm $dir $N $tag $x ${x1[$SGE_TASK_ID-1]} $nbins_t $max_eig $n_threads 


# dir=test_run

# n_threads=50
# max_eig=225
# nbins_t=75

# N=5
# tag=1
# x=(0.0 0.0 0.0 0.0  0.001 0.001 0.001 0.001)
# x1=(0.0 0.05 0.1 0.15  0.0 0.05 0.1 0.15)

# ./nm $dir $N $tag $x ${x1[$SGE_TASK_ID-1]} $nbins_t $max_eig $n_threads 



# ./nm test_run 5 1 0.001 0.1 100 144 50 


# dir=how_Vkl_scales_with_N,threads=100

# N=(5 7 9 11  5 7 9 11  5 7 9 11  5 7 9 11  5 7 9 11)
# tag=(3 4 5 6  3 4 5 6  3 4 5 6  3 4 5 6  3 4 5 6)
# max_eig=(144 144 144 144  196 196 196 196  225 225 225 225  289 289 289 289  400 400 400 400)
# n_threads=100
# nbins_t=1
# x=0.5
# x1=0.6

# ./nm $dir ${N[$SGE_TASK_ID-1]} ${tag[$SGE_TASK_ID-1]} $x $x1 $nbins_t ${max_eig[$SGE_TASK_ID-1]} $n_threads 