#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/ # root directory for project (non-scratch)
SNAKE=/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/analysis/cluster_inbreeding
NODES=1028 # max number of cluster nodes
WAIT=600 # lag for system

snakemake \
	--snakefile $SNAKE/run_snake_params.py \
	--configfile $SNAKE/config.yaml \
	--directory $ROOT \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--cluster $SNAKE/launch.py \
	-j $NODES \
#	--dryrun -p
