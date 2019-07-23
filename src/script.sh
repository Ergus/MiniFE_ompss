#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

module purge
module load bsc/1.0
module load gcc/7.2.0

module load hwloc/1.11.8
module load intel/2018.1 impi/2018.1 mkl/2018.1
module load EXTRAE/3.6.1
module load boost/1.64.0-mpi
module load nanos6/cluster_impi
module load mcxx/git_impi

ulimit -c unlimited

# Environment variables
export NANOS6_CPU_SCHEDULER=fifo
export NANOS6_COMMUNICATION=mpi-2sided
export NANOS6_DISTRIBUTED_MEMORY=128G
export NANOS6_LOCAL_MEMORY=10G

export NANOS6_SCHEDULER=cluster-locality
export NANOS6_COMMUNICATION=mpi-2sided
export NANOS6=optimized

nx=$1
np=$2
numboxes=$3

echo -e "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo -e "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo -e "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo -e "# QOS: ${SLURM_JOB_QOS}"
echo -e "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
env | grep NANOS6 | sed -e 's/^#*/# /'
echo -e "# ======================================\n"

COMMAND="mpirun -np ${np} ./miniFE.x nx=${nx} numboxes=${numboxes}"
echo -e "# Command: " ${COMMAND}
${COMMAND}
