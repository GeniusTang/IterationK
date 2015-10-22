#!/bin/sh

while getopts k:t:m: opts; do
    case $opts in
        k) K=$OPTARG ;;
        t) T=$OPTARG ;;
        m) M=$OPTARG
    esac
done


cd /home/rcf-40/kujin/panasas/IterationK

rm -rf /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_k${K}_${M}_rf
touch /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_k${K}_${M}_rf
rm -rf /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_iter_k${K}_${M}_rf
touch /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_iter_k${K}_${M}_rf
for j in {1..100};do
    cd /home/rcf-40/kujin/panasas/IterationK/k$K/T$T 
    FILEPATH=/home/rcf-40/kujin/panasas/Examples/T$T/rep$j/
    python /home/rcf-40/kujin/panasas/IterationK/iterationk.py -i ${FILEPATH}file -k $K -t L -m ${M} -o ${FILEPATH}${M}_k${K}_iter_ --iter
    python /home/rcf-40/kujin/panasas/IterationK/iterationk.py -i ${FILEPATH}file -k $K -t L -m ${M} -o ${FILEPATH}${M}_k${K}_ 
    echo "(((A,B),(C,D)),((G,H),(E,F)));" >> ${FILEPATH}${M}_k${K}_iter_outtree
    echo "(((A,B),(C,D)),((G,H),(E,F)));" >> ${FILEPATH}${M}_k${K}_outtree
    echo -e "${FILEPATH}${M}_k${K}_iter_outtree\nR\nD\n1\n2\nA\nS\nY" > input1
    echo -e "${FILEPATH}${M}_k${K}_outtree\nR\nD\n1\n2\nA\nS\nY" > input2
    rm -rf outfile
    treedist < input2 > screenout 
    cat outfile >> /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_k${K}_${M}_rf
    rm -rf outfile 
    treedist < input1 > screenout 
    cat outfile >> /home/rcf-40/kujin/panasas/IterationK/RF/T${T}_iter_k${K}_${M}_rf
done
