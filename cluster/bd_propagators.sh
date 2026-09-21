#!/bin/bash
#$ -S /bin/bash
#$ -q *
#$ -N g2p5-50
#$ -l h_rt=50:00:00
#$ -t 1-6:1
#$ -cwd

output_dir="g2p2"

type=(m m nm nm nm nm)
confinement=(1.0 1.0 1.0 1.0 0.05 0.05)
# M=10000    # ten thousand
# M=25000    # 25 thousand
M=50000    # fifty thousand

N=5
tag=(1 3 1 3 1 3)

# N=11
# tag=(1 6)

# N=15
# tag=(1 8)

# N=21
# tag=(1 11)

# N=31
# tag=(1 16)

# N=51
# tag=(1 26)

~/bin/julia g2p.jl $output_dir ${type[$SGE_TASK_ID-1]} $N ${tag[$SGE_TASK_ID-1]} ${confinement[$SGE_TASK_ID-1]} $M



