#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

NODES=(64 128 256 512 1024 2048 4096 8192 16384)

CLIQUES=(8 12 16 22 32 45 64 90 128)

for N in "${NODES[@]}"; 
do
    
    ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 5)) --seed 3461 --weight 10 --randomize > ../test/graphs/weighted/erdos_renyi_Ntimes5/er_${N}.wel
    # ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 60)) --weight 1 > ../test/graphs/unweighted/erdos_renyi_Ntimes60/er_${N}.wel
    # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 8 --weight 1 > ../test/graphs/unweighted/watts_strogatz_8_0.3/ws_${N}.wel
    # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 1 > ../test/graphs/unweighted/watts_strogatz_24_0.3/ws_${N}.wel
    # ${DIR}/generate.py 'barabasi_albert_graph(N,3)' ${N} --weight 1 > ../test/graphs/unweighted/barabasi_albert_3/ba_${N}.wel
    # ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 1 > ../test/graphs/unweighted/barabasi_albert_9/ba_${N}.wel

    # ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 15)) --weight 10 --randomize > ../test/graphs/weighted/erdos_renyi_Ntimes15/er_${N}.wel
    # ${DIR}/generate.py 'dense_gnm_random_graph(N,M)' ${N} -m $((N * 60)) --weight 10 --randomize> ../test/graphs/weighted/erdos_renyi_Ntimes60/er_${N}.wel
    # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 8 --weight 10 --randomize > ../test/graphs/weighted/watts_strogatz_8_0.3/ws_${N}.wel
    # ${DIR}/generate.py 'connected_watts_strogatz_graph(N,K,0.3)' ${N} -k 24 --weight 10 --randomize > ../test/graphs/weighted/watts_strogatz_24_0.3/ws_${N}.wel
    # ${DIR}/generate.py 'barabasi_albert_graph(N,3)' ${N} --weight 10 --randomize > ../test/graphs/weighted/barabasi_albert_3/ba_${N}.wel
    # ${DIR}/generate.py 'barabasi_albert_graph(N,9)' ${N} --weight 10 --randomize > ../test/graphs/weighted/barabasi_albert_9/ba_${N}.wel
done

# for C in "${CLIQUES[@]}";
# do
#     ${DIR}/generate.py 'ring_of_cliques(NC,NC)' 1 -nc ${C} --weight 1 > unweighted/roc_$((C*C)).wel
#     ${DIR}/generate.py 'ring_of_cliques(NC,NC)' 1 -nc ${C} --weight 10 --randomize > weighted/roc_$((C*C)).wel

# done