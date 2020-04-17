#! /usr/bin/env python3
"""
Master script for running MVN RF
"""

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
clst_inb_paramDIR = config["paramdir"]
clst_inb_outDIR = config["outdir"]

final_target = []
for i in parampath:
	final_target.append( os.path.join(clst_inb_outDIR, "{}.RDS".format(i)) )


rule all:
	input: final_target

rule params_out:
	input: clst_inb_paramDIR + "{params}.RDS"
	output: clst_inb_outDIR + "{params}.RDS",
	log: clst_inb_outDIR + "{params}_log.Rout",
	shell:
		r"""
		Rscript --max-ppsize=500000 --vanilla \
			R/clst_inb_coeff_run.R \
			--mastermap {input} \
			-O {output} \
			>& {log}
		"""
