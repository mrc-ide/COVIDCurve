#! /bin/bash
#SNAKE=/proj/ideel/meshnick/users/NickB/Projects/CurveAware/ignore/deploy/cluster_deploy
SNAKE=/Users/nickbrazeau/Documents/GitHub/CurveAware/ignore/deploy/cluster_deploy
NODES=1028 # max number of cluster nodes
WAIT=60 # lag for system

snakemake \
	--snakefile $SNAKE/run_snake_mcmc.py \
	--configfile $SNAKE/config.yaml \
	--directory $SNAKE \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--cluster $SNAKE/launch.py \
	-j $NODES \
	--dryrun -p
