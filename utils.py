

def load_rep_genome_cluster(gtdb_cluster_loc):
	"""
	returns dictionary of 
	taxonomy : [representative genome, [all genomes]]

	requires sp_clusters.tsv from 
	https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/
	"""
	gtdb_genomes = {}

	with open(gtdb_cluster_loc, "r") as f:
		next(f)
		for row in f:
			row = row.replace("\n","").split("\t")

			rep_genome = row[0]
			taxonomy = row[1].split("__")[1].lower()
			all_genomes = row[9].split(",")

			gtdb_genomes[taxonomy] = [rep_genome, all_genomes]
	f.close()

	return gtdb_genomes


def get_species_from_rankofinterest(tree, rank_of_interest):
	"""
	"""
	taxids_from_rank = []
	for taxonomy in tree.nodes:
		if tree[taxonomy].rank == rank_of_interest:
			taxids_from_rank.append(taxonomy)

	taxid_rank_dict = {}
	for rank_taxonomy in taxids_from_rank:
		taxid_rank_dict[rank_taxonomy] = []

		child_taxids = tree.descend(rank_taxonomy)
		for c in child_taxids:
			if tree[c].rank == "species":
				taxid_rank_dict[rank_taxonomy].append(c)

	return taxid_rank_dict