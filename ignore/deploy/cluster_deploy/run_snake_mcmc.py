#! /usr/bin/env python3

from __future__ import print_function

import os
import sys
import yaml
import re
from collections import defaultdict

def load_run_metadata(f):
	""" Get params from tab-separated file."""
	params = list()
	with open(f) as rfile:
		for line in rfile:
			if line.startswith("#"):
				continue
			line = line.strip().split("\t")
			params.append(line[0])
	return params


## read run manifest and neccessary supports
parampath = load_run_metadata(config["manifest"])
paramDIR = config["paramdir"]
outDIR = config["outdir"]

final_target = []
for i in parampath:
	final_target.append( os.path.join(outDIR, "{}.rds".format(i)) )


rule all:
	input: final_target

rule params_out:
	input: paramDIR + "{params}.rds"
	output: outDIR + "{params}.rds",
	log: outDIR + "{params}_log.Rout",
	shell:
		r"""
		Rscript --max-ppsize=500000 --vanilla \
			wrap_mcmc_cluster_runs.R \
			--mastermap {input} \
			-O {output} \
			>& {log}
		"""
