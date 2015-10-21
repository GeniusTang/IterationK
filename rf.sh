#!/bin/bash
QSUB="qsub -q cmb -l nodes=1:ppn=4 -l mem=20000m -l walltime=100:00:00"

cd /home/rcf-40/kujin/panasas/IterationK


for i in {3..4};do
    echo "./rfcompute.sh -k $i" | $QSUB
done

