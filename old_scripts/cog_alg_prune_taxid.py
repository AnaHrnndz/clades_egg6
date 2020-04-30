from ete3  import SeqGroup

## I need to remove taxid that are not leaves in guide tree from the alg to run raxml (to optimize branch lengths in the bifurcating tree)
taxids = []
with open('/home/plaza/projects/eggnog6/guide_tree_ncbi/final_taxid_list') as final_taxid:
    for line in final_taxid:
        if line.strip():
            sp = str(line.strip())
            taxids.append(sp)
print('len taxid', len(taxids))

alg = SeqGroup('/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/COGs_seqs.fasta.final_tree.fa')
alg_clean= open('/home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/prune_COGs_seqs.fasta.final_tree.fa', 'w')

for num, (name, seq, _) in enumerate(alg):
    if str(name) in taxids:
        print('>%s\n%s' %(name, seq), file = alg_clean)

alg_clean.close()



