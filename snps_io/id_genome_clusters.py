from __future__ import division

import sys, os
import argparse, operator

import numpy as np

from time import time

class GenomeCluster:
	def __init__(self, max_d):
		self.max_d  = max_d

		self.genomes = dict()
		self.links = dict()

		self.tag_genome = None

	def size(self):
		return len(genomes)

	def add(self, genome1, genome2, d):
		if genome1 not in self.genomes:
			self.genomes[genome1] = 1
		else:
			self.genomes[genome1] = self.genomes[genome1] + 1

		if genome2 not in self.genomes:
			self.genomes[genome2] = 1
		else:
			self.genomes[genome2] = self.genomes[genome2] + 1

		link = "{}|{}".format(genome1, genome2)
		if genome1 > genome2:
			link = "{}|{}".format(genome2, genome1)

		#sys.stderr.write("{} {}\n".format(genome1, genome2))

		if link not in self.links:
			self.links[link] = d
		else:
			#assert True
			assert self.links[link] == d
			#sys.stderr.write("{} {}\n".format(genome1, genome2))
			

	def merge(self, cluster2, genome1, genome2, d):
		if self.contains(genome1) and cluster2.contains(genome2):
			self.add(genome1, genome2, d)
			for link in cluster2.links.keys():
				genomes = link.split("|")
				self.add(genomes[0], genomes[1], cluster2.links[link])
		else:
			sys.exit("{} is not in cluster1 or {} is not cluster2, nullify the basis for merging".format(genome1, genome2))

	def is_empty(self):
		return len(self.genomes.keys()) == 0

	def contains(self, genome):
		return genome in self.genomes

	def id_tag_genome(self):
		# no snps
		if len(self.genomes) == 0:
			sys.exit("\nError: no genomes on cluster: cannot id tag genome\n")
		# one snp
		elif len(self.genomes) == 1:
			self.tag_genome = self.genomes.keys()[0]
		else:
			tmp_min = 0
			tmp_genome = None
			for genome in self.genomes.keys():
				if self.genomes[genome] > tmp_min:
					tmp_min = self.genomes[genome]
					tmp_genome = genome
			self.tag_genome = tmp_genome

		return self.tag_genome

	def fmtout(self):
		sorted_tuples = sorted(self.genomes.items(), key=operator.itemgetter(1), reverse=True)
		sorted_genomes = [genome_tuple[0] for genome_tuple in sorted_tuples]

		return "* {} {}".format(self.tag_genome, " ".join(sorted_genomes))

	def fmtout_all(self):
		fmt_str = "{}\n".format(self.fmtout())

		for link in self.links.keys():
			fmt_str += "- {} {}\n".format(link, self.links[link])

		return fmt_str

def search_genome_clusters(dist_path, max_d):
	sys.stderr.write("[clustering] start\n")

	genome_clusters = []
	genome_lookup = dict()

	with open(dist_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split('\t')
			genome1, genome2, d = items[0], items[1], float(items[2])

			if genome1 >= genome2 or d > max_d:
				#sys.stderr.write("{} {}\n".format(genome1, genome2))
				continue
			# sys.stderr.write("{} {}\n".format(genome1, genome2))

			if genome1 not in genome_lookup and genome2 not in genome_lookup:
				new_cluster = GenomeCluster(max_d)
				new_cluster.add(genome1, genome2, d)
				genome_clusters.append(new_cluster)
				genome_lookup[genome1] = len(genome_clusters) - 1
				genome_lookup[genome2] = len(genome_clusters) - 1
			elif genome1 in genome_lookup and genome2 not in genome_lookup:
				cluster_indx = genome_lookup[genome1]
				genome_lookup[genome2] = cluster_indx
				genome_clusters[cluster_indx].add(genome1, genome2, d)
			elif genome1 not in genome_lookup and genome2 in genome_lookup:
				cluster_indx = genome_lookup[genome2]
				genome_lookup[genome1] = cluster_indx
				genome_clusters[cluster_indx].add(genome1, genome2, d)
			else:
				if genome_lookup[genome1] == genome_lookup[genome2]:
					pass
				else:
					cluster_indx1 = genome_lookup[genome1]
					cluster_indx2 = genome_lookup[genome2]
					genome_clusters[cluster_indx1].merge(genome_clusters[cluster_indx2], genome1, genome2, d)

					for genome in genome_clusters[cluster_indx2].genomes:
						genome_lookup[genome] = cluster_indx1

					genome_clusters[cluster_indx2] = None

	sys.stderr.write("[clustering] done\n")
	sys.stderr.write("[clustering] {} genomes have been included in clusters\n".format(len(genome_lookup.keys())))

	good_clusters = verify_clusters(genome_clusters, genome_lookup)

	return good_clusters, len(genome_lookup.keys())

def verify_clusters(genome_clusters, genome_lookup):
	for genome in genome_lookup.keys():
		assert genome in genome_clusters[genome_lookup[genome]].genomes

	good_clusters = []
	for i, cluster in enumerate(genome_clusters):
		if cluster is not None:
			cluster.id_tag_genome()
			good_clusters.append(cluster)

			for genome in cluster.genomes:
				assert genome_lookup[genome] == i

	return good_clusters

def output_clusters(good_clusters, output_path="/dev/stdout"):
	if output_path is not None:
		with open(output_path, 'w') as fh:
			for gcluster in good_clusters:
				fh.write(gcluster.fmtout_all())

def build_genome_blocks(dist_path, total_n, critical_n=100, max_d=0.01, end_d=0.000001, range_factor=1.2, output_path=None):
	optimal_d = 0
	optimal_n = 0
	optimal_clusters = []

	upper_cap = critical_n * range_factor

	genome_clusters, clust_n = search_genome_clusters(dist_path, max_d)
	
	tag_n = total_n - clust_n + len(genome_clusters)

	firstcut_exit = False
	if tag_n > upper_cap:
		print("Program will continue with a non-optimal number ({}) of genomes. Perhaps try a higher cutoff (current {})".format(str(tag_n), str(max_d)))
		optimal_d = max_d
		optimal_n = tag_n
		optimal_clusters = genome_clusters
		firstcut_exit = True
	elif tag_n >= critical_n and tag_n <= upper_cap:
		# perfect scenario on exit
		optimal_d = max_d
		optimal_n = tag_n
		optimal_clusters = genome_clusters
	else:
		# determine lower bound
		min_d = max_d

		print("[Searching lower cap]")
		while min_d >= end_d and tag_n < critical_n:
			min_d = min_d / 10

			genome_clusters, clust_n = search_genome_clusters(dist_path, min_d)
			tag_n = total_n - clust_n + len(genome_clusters)

			print("\t{}: {} tag genomes".format(min_d, tag_n))

		print("[End earching]")

		# binary search into critical range
		print("[Searching optimal d-cut]")
		if min_d < end_d and tag_n < critical_n:
			print("Program cannot reach the number ({}) of genomes required for core-genome SNP calling.")
			print("Proceeding with orginal set of genomes")

			optimal_d = None 
			optimal_n = None 
			optimal_clusters = None 
		else:
			left_d = max_d
			right_d = min_d
			mid_point = int((upper_cap + critical_n) / 2)

			delta_d = 1 # arbitary value; does not matter

			while delta_d > 0.0000001 and (tag_n > upper_cap or tag_n < critical_n):
				cur_d = (left_d + right_d) / 2

				genome_clusters, clust_n = search_genome_clusters(dist_path, cur_d)
				tag_n = total_n - clust_n + len(genome_clusters)

				if tag_n > mid_point:
					right_d = cur_d
				else:
					left_d = cur_d

				delta_d = abs(left_d - right_d)
					
				print("\tsearching space [ {} , {} ]".format(left_d, right_d))
				print("\tcurrent d-cut: {}".format(cur_d))
				print("\tcurrent no of tags: {}".format(tag_n))

				if tag_n >= critical_n:
					delta_1 = abs(tag_n - mid_point)
					delta_2 = abs(optimal_n - mid_point)

					if delta_1 < delta_2:
						optimal_d = cur_d
						optimal_n = tag_n
						optimal_clusters = genome_clusters
					else:
						pass
				else:
					pass


		print("[Searching optimal d-cut]")

		if tag_n < critical_n and optimal_n < critical_n:
			print("Program cannot reach the number ({}) of genomes required for core-genome SNP calling.")
			print("Proceeding with orginal set of genomes. Or try higher MAF")

			optimal_d = None 
			optimal_n = None 
			optimal_clusters = None 

	return optimal_clusters, optimal_d, optimal_n, firstcut_exit 
