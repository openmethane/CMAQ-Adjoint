#!/bin/bash
#PBS -P q90
#PBS -q normal
#PBS -l walltime=20:00,mem=16GB
#PBS -l ncpus=16
#PBS -l wd

source /home/563/spt563/mods/module_cmaq.sh

#cmaq_adj crashes when environment variable set with no value.
unset PBS_NCI_IMAGE

./run.adj.fwd.bnmk &> fwd.bnmk.log
./run.adj.bwd.bnmk &> bwd.bnmk.log
