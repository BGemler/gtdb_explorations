from config import *
from utils import load_rep_genome_cluster, get_species_from_rankofinterest
from outwriter import write_out_taxid_splice


def main():
	"""
	"""
	# get representative genome dict for each taxonomy
	gtdb_genomes = load_rep_genome_cluster(gtdb_cluster_loc)

	# get dictionary of GTDB phylums, with all species associated with each phylum
	taxid_phylum_dict = get_species_from_rankofinterest(tree, splice_rank_of_interest)

	# write out CSV of each phylum with stats of interest
	write_out_taxid_splice(tree, \
														splice_rank_of_interest, \
														taxid_phylum_dict, gtdb_genomes)


	return


main()