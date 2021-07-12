import csv


def convert_list_to_out(list_in):
	str_out = ""
	for l in list_in:
		str_out = str_out + str(l) + "\n"
	str_out = str_out[:-1]

	return str_out


def write_out_taxid_splice(tree, \
														splice_rank_of_interest, \
														taxid_rank_dict, gtdb_genomes):
	"""
	"""
	splice_write_out_loc = "results/" + splice_rank_of_interest + "-gtdb_slice.csv"
	with open(splice_write_out_loc, "w") as f:
		out = csv.writer(f)
		out.writerow([
			"Domain Taxonomy", "Taxonomy Splice at Rank:" + splice_rank_of_interest, \
					"Number of Species in Splice", \
					"Total Number of Genomes from Species in Splice", \
					"Top Species by Number of Genomes"
			])

		for taxid in taxid_rank_dict:
			domain_taxonomy = tree.ascend(taxid)[-1]

			species = taxid_rank_dict[taxid]
			num_species = len(species)

			tot_num_genomes = 0
			rep_species_genomes = []
			species_count_list = []
			for s in species:
				rep_genome, tot_genomes = gtdb_genomes[s]

				tot_num_genomes += len(tot_genomes)
				rep_species_genomes.append(rep_genome) 

				species_count_list.append([s, len(tot_genomes)])
			rep_species_genomes = rep_species_genomes[:-1]

			species_count_list = sorted(species_count_list, key = lambda x:x[1], reverse = True)
			down_species_count_list = species_count_list[:5] if len(species_count_list) > 5 else species_count_list

			out.writerow([domain_taxonomy, taxid, num_species, tot_num_genomes, convert_list_to_out(down_species_count_list)])


	return