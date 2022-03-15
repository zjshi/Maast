from __future__ import division

import sys, os
import argparse, operator

import numpy as np

from time import time

def read_tags(tag_paths):
	tag_map = dict()

	for tag_genome in tag_paths:
		tag_map[tag_genome] = 0
	
	return tag_map

def calc_tag_weights(tag_map, dist_path):
	sys.stderr.write("[clustering] start\n")

	with open(dist_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split('\t')
			genome1, genome2, d = items[0], items[1], float(items[2])

			if genome1 >= genome2:
				#sys.stderr.write("{} {}\n".format(genome1, genome2))
				continue
			# sys.stderr.write("{} {}\n".format(genome1, genome2))

			if genome1 in tag_map and genome2 in tag_map:
				tag_map[genome1] += d
				tag_map[genome2] += d

	sys.stderr.write("[clustering] done\n")

	return tag_map

def centroid_from_map(tag_map):
	centroid = None

	for tag in tag_map.keys():
		if centroid is None:
			centroid = tag
		else:
			if tag_map[tag] < tag_map[centroid]:
				centroid = tag

	return centroid

def identify(tag_paths, dist_path):
	tag_map = read_tags(tag_paths)
	tag_map = calc_tag_weights(tag_map, dist_path)

	centroid = centroid_from_map(tag_map)

	return centroid
