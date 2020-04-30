from ete3 import Tree
#We have to remove 1458427 and 595593 taxid from concatenate tree 
#because we dont have those taxid in reference tree build with ncbi taxa

in_tree = Tree('/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/COGs_seqs.fasta.final_tree.nw')
#out_tree = Tree('/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/prune_COGs_seqs.fasta.final_tree.nw')

print(len(in_tree))
n1 = in_tree.search_nodes(name="1458427")[0]
n2 = in_tree.search_nodes(name="595593")[0]

print (n1)
print (n2)

removed_node = n1.detach()
removed_node = n2.detach()

print(len(in_tree))
in_tree.resolve_polytomy()

in_tree.write(format=1, outfile='/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/bifurcating_prune_COGs_seqs.fasta.final_tree.nw')