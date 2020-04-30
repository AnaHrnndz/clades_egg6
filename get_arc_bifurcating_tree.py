from ete3 import Tree


#Clean the tree to contain only taxids in leaf nodes
ref_tree = Tree('/home/plaza/projects/eggnog6/ncbi_tree/Archaea_ncbi.nw')

## randomly resolve politomies in the reference tree
ref_tree.resolve_polytomy()


final_taxid_list = []
for n in ref_tree.traverse():
    if n.is_leaf():
        n.name = n.name
        final_taxid_list.append(n.name)
    else:
        n.name = ''

final_taxid=open('/home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/arc_final_taxid_list','w')
final_taxid.write('\n'.join(final_taxid_list))
final_taxid.close()

ref_tree.write(outfile='/home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/arc_bif_ncbi_tree.nw', format = 1)
