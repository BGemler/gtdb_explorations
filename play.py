import gtdb_tree
print("loading gtdb taxtree")
tree = gtdb_tree.TaxonomyTree('resources/gtdb_taxdmp')
print("finished loading gtdb taxtree\n\n")



taxid_of_interest = "coriobacteriia"

print(taxid_of_interest)
print(tree[taxid_of_interest])

print("\n\nascending:")
ascension = tree.ascend(taxid_of_interest)
for a in ascension:
	print(a)


print("\n\ndescending:")
descension = tree.descend(taxid_of_interest)
for d in descension:
	print(d)


print("\n\n\n\n\nbacteria:")

print(tree["bacteria"])
