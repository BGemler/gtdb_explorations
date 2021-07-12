import gtdb_tree
print("loading gtdb taxtree")
tree = gtdb_tree.TaxonomyTree('resources/gtdb_taxdmp')
print("finished loading gtdb taxtree\n\n")

# paths
gtdb_cluster_loc = "resources/sp_clusters.tsv"


# variables
splice_rank_of_interest = "phylum"