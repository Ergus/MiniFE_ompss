#!/bin/bash

np_array=(1 2 4 8)
nx_array=(64 128 256 512 1024)
numboxes_array=(4 8 16 32 64)

now=$(date +%F_%H-%M-%S)
resdir="results/minife_${now}"

mkdir -p ${resdir}

echo "Directory: ${resdir}"
echo "Submitting: group"
echo "nodes: ${np_array[*]}"
echo "numboxes: ${numboxes_array[*]}"
echo "dim: ${nx_array[*]}"

for nx in ${nx_array[@]}; do
    for np in ${np_array[@]}; do
	for numboxes in ${numboxes_array[@]}; do

	    printf "Submitting combination nx: %s np: %s numboxes: %s\n" \
		   $nx $np $numboxes
	    jobname="mfe_${nx}_${np}_${numboxes}"
	    filename="${resdir}/${jobname}"

 	    sbatch --ntasks=${np} \
		   --qos="debug" \
 		   --job-name=${jobname} \
 		   --output="${resdir}/%x_%2a_%j.out" \
 		   --error="${resdir}/%x_%2a_%j.err" \
 		   ./script.sh ${nx} ${np} ${numboxes}
	done
    done
done
