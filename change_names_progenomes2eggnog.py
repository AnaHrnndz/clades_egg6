from ete3 import SeqGroup


alg = SeqGroup('/home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree/cog_all-alg_concat_default-fasttree_default/Bact_COGs_seqs.faa.final_tree.fa')
alg_clean = open('/home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree/cog_all-alg_concat_default-fasttree_default/clean_Bact_COGs_seqs.faa.final_tree.fa', 'w')

eggnog62progenomes2 = {}
progenomes2eggnog = {}
with open('/home/plaza/projects/eggnog6/raw_data/eggnogv6_2_progenomesv2.tsv') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            info = line.split()
            taxid_egg6 = info[0].strip()
            taxid_pro2 = info[1].strip()
            eggnog62progenomes2[int(taxid_egg6)] = int(taxid_pro2)
            progenomes2eggnog[int(taxid_pro2)] = int(taxid_egg6)


for num, (name, seq, _) in enumerate(alg):
    name=int(name)
    taxid = progenomes2eggnog[name]
    alg_clean.write('>%s\n%s\n' %(taxid, seq))

alg_clean.close()

