#----------------------------------------------------------------------------
# cNMF_exec1.sh
#----------------------------------------------------------------------------
#!/bin/sh
i=$1

cnmf prepare \
--output-dir ./output/cNMF/result1/ \
--name $i \
-c ./output/cNMF/input_data1/counts_${i}.txt \
-k 4 5 6 7 8 9 \
--n-iter 50 --seed 14

cnmf factorize \
--output-dir ./output/cNMF/result1/ \
--name $i \
--worker-index 0 --total-workers 1

cnmf combine \
--output-dir ./output/cNMF/result1/ \
--name $i

cnmf k_selection_plot \
--output-dir ./output/cNMF/result1/ \
--name $i

for n in {4..9}; do
cnmf consensus \
--output-dir ./output/cNMF/result1/ \
--name $i \
--components $n --local-density-threshold 0.2 --show-clustering
done
