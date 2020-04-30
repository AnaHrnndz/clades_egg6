from ete3 import Tree


#Clean the tree to contain only taxids in leaf nodes
ref_tree = Tree('/home/plaza/projects/eggnog6/guide_tree_ncbi/freeze12_representative_sp_tree_oldDB.nw')

## randomly resolve politomies in the reference tree
ref_tree.resolve_polytomy()


final_taxid_list = []
for n in ref_tree.traverse():
    if n.is_leaf():
        n.name = n.name
        final_taxid_list.append(n.name)
    else:
        n.name = ''

final_taxid=open('/home/plaza/projects/eggnog6/guide_tree_ncbi/final_taxid_list','w')
final_taxid.write('\n'.join(final_taxid_list))
final_taxid.close()

ref_tree.write(outfile='/home/plaza/projects/eggnog6/guide_tree_ncbi/post_process_freeze12_representative_sp_tree_oldDB.nw', format = 1)
