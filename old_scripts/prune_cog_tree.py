from ete3 import Tree

#prune COG tree for keeping only taxid that are leaves in guide tree (ncbi taxonomy)
taxids = []
with open('/home/plaza/projects/eggnog6/guide_tree_ncbi/final_taxid_list') as final_taxid:
    for line in final_taxid:
        if line.strip():
            sp = str(line.strip())
            taxids.append(sp)
print('len taxid', len(taxids))

in_tree = Tree('/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/COGs_seqs.fasta.final_tree.nw')
#in_tree.resolve_polytomy()

print (len(in_tree))
#print ('resolved politomy')

in_tree.prune(taxids)
print ('prune tree')
print (len(in_tree))
in_tree.write(format=1, outfile='/home/plaza/projects/eggnog6/prune_cog_tree/bifurcating_prune_COGs_seqs.fasta.final_tree.nw')