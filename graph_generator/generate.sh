#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

NODES=(64 128 256 512 1024 2048 4096 8192 16384)

CLIQUES=(8 12 16 22 32 45 64 90 128)

PARTITION_SIZES=(21 42 85 170 341 682 1365 2730 5461)
# PARTITION_SIZES=(5461)

# P_OUT=(0.00396825, 0.00198413, 0.000980392, 0.000490196, 0.000244379, 0.00012219, 0.0000610501, 0.000030525, 0.0000152597)

# for N in "${NODES[@]}"; 
# do
    
#     # ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --seed 3461 --weight 10 --randomize > ../test/graphs/weighted/erdos_renyi_Ntimes5/er_${N}.wel
#     # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 1 > ../test/graphs/unweighted/watts_strogatz_24_0.3/ws_${N}.wel
#     # ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 1 > ../test/graphs/unweighted/barabasi_albert_9/ba_${N}.wel

#     # ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --weight 10 --randomize > ../test/graphs/weighted/erdos_renyi_Ntimes15/er_${N}.wel
#     # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 10 --randomize > ../test/graphs/weighted/watts_strogatz_24_0.3/ws_${N}.wel
#     # ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 10 --randomize > ../test/graphs/weighted/barabasi_albert_9/ba_${N}.wel
# done

# for C in "${CLIQUES[@]}";
# do
#     ${DIR}/generate.py 'ring_of_cliques(NC,NC)' 1 -nc ${C} --weight 1 > unweighted/roc_$((C*C)).wel
#     ${DIR}/generate.py 'ring_of_cliques(NC,NC)' 1 -nc ${C} --weight 10 --randomize > weighted/roc_$((C*C)).wel

# done

for ((i = 0 ; i < ${#PARTITION_SIZES[@]} ; i++));
do
    TMP=$((PARTITION_SIZES[i] * 3))
    TMP1=`python3 -c "print(1/(4*${TMP}))"`
    ${DIR}/generate.py 'random_partition_graph(PL, PIN, POUT)' 1 -pl ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} -pin 0.25 -pout ${TMP1} --weight 1 > ../test/graphs/unweighted/random_partition_graph_0.25_1div4n/rp_${TMP}.wel

done