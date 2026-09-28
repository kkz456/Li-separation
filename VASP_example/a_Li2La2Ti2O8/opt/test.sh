#!/bin/bash
#SBATCH -J xxx
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p wzacnormal
module purge
module load vasp-6.4.2-intelmpi2017_ioptcell

export MKL_DEBUG_CPU_TYPE=5 #���ٴ���
export MKL_CBWR=AVX2 #ʹcpuĬ��֧��avx2
export I_MPI_PIN_DOMAIN=numa #�ڴ�λ����cpuλ�ð󶨣������ڴ��ȡ�������ڴ����Ҫ��ߵļ�����������

srun --mpi=pmi2 vasp_std
