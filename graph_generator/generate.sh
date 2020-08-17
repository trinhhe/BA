#!/usr/bin/env bash
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=/cluster/home/trinhhe/gapbs/graph_generator

NODES=(512 1024 2048 4096 8192 12288 16384 20480 24576 )
# NODES=(28672 33280)
# NODES=(64 128)

# PARTITION_SIZES=(500 1000 1500)
# PARTITION_SIZES=(5461)

NUMBER_OF_GROUPS=(4 8 16 32 64 96 128 160 192 )
# NUMBER_OF_GROUPS=(224 260)
# P_OUT=(0.00396825, 0.00198413, 0.000980392, 0.000490196, 0.000244379, 0.00012219, 0.0000610501, 0.000030525, 0.0000152597)

for N in "${NODES[@]}"; 
do
    
    ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --weight 1 > ../../graphs/unweighted/erdos_renyi_Ntimes5/er_${N}.graph
    ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 1 > ../../graphs/unweighted/watts_strogatz_24_0.3/ws_${N}.graph
    ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 1 > ../../graphs/unweighted/barabasi_albert_9/ba_${N}.graph

    ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --weight 100 --sigma 5 --randomize_normal > ../../graphs/weighted_100_5/erdos_renyi_Ntimes5/er_${N}.graph
    ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 100 --sigma 5 --randomize_normal > ../../graphs/weighted_100_5/watts_strogatz_24_0.3/ws_${N}.graph
    ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 100 --sigma 5 --randomize_normal > ../../graphs/weighted_100_5/barabasi_albert_9/ba_${N}.graph

    ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --weight 100 --sigma 30 --randomize_normal > ../../graphs/weighted_100_30/erdos_renyi_Ntimes5/er_${N}.graph
    ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 100 --sigma 30 --randomize_normal > ../../graphs/weighted_100_30/watts_strogatz_24_0.3/ws_${N}.graph
    ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 100 --sigma 30 --randomize_normal > ../../graphs/weighted_100_30/barabasi_albert_9/ba_${N}.graph
done


# for ((i = 0 ; i < ${#PARTITION_SIZES[@]} ; i++));
# do
#     TMP=$((PARTITION_SIZES[i]*4))
#     TMP1=`python3 -c "print(1/(${PARTITION_SIZES[i]}))"`
#     TMP2=`python3 -c "print(1/(4*${TMP}))"`
#     ${DIR}/generate.py 'random_partition_graph(PL, PIN, POUT)' 1 -pl ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} -pin ${TMP1} -pout ${TMP2} --weight 1 > ../../graphs/unweighted/random_partition_graph_1divn_1div4n/rp_${TMP}.graph
#     ${DIR}/generate.py 'random_partition_graph(PL, PIN, POUT)' 1 -pl ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} ${PARTITION_SIZES[i]} -pin ${TMP1} -pout ${TMP2} --weight 100 --randomize_normal > ../../graphs/weighted_100_5/random_partition_graph_1divn_1div4n/rp_${TMP}.graph

# done

for ((i = 0 ; i < ${#NUMBER_OF_GROUPS[@]} ; i++));
do
    TMP=$((NUMBER_OF_GROUPS[i]*128))
    TMP2=`python3 -c "print(1/(4*${TMP}))"`
    ${DIR}/generate.py 'planted_partition_graph(L, K, PIN, POUT)' 1  -l ${NUMBER_OF_GROUPS[i]} -k 128 -pin 0.5 -pout ${TMP2} --weight 1 > ../../graphs/unweighted/planted_partition_128_0.5_1div4n/pp_${TMP}.graph
    ${DIR}/generate.py 'planted_partition_graph(L, K, PIN, POUT)' 1  -l ${NUMBER_OF_GROUPS[i]} -k 128 -pin 0.5 -pout ${TMP2} --weight 100 --sigma 25 --randomize_normal > ../../graphs/weighted_100_25/planted_partition_128_0.5_1div4n/pp_${TMP}.graph
    ${DIR}/generate.py 'planted_partition_graph(L, K, PIN, POUT)' 1  -l ${NUMBER_OF_GROUPS[i]} -k 128 -pin 0.5 -pout ${TMP2} --weight 100 --sigma 5 --randomize_normal > ../../graphs/weighted_100_5/planted_partition_128_0.5_1div4n/pp_${TMP}.graph

done