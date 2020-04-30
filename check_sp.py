from ete3 import NCBITaxa

ncbi = NCBITaxa('/home/plaza/projects/eggnog6/ncbi_taxa_db/taxa.sqlite')


ori_arc_sp = []
with open('/home/plaza/projects/eggnog6/ncbi_tree/arch_sp_tree.tsv') as ori_arc:
    for line in ori_arc:
        if line.strip() and not line.startswith('#'):
            sp = line.strip()
            ori_arc_sp.append(sp)

print (len(ori_arc_sp))


mems = []
with open('/home/plaza/projects/eggnog6/arc_newXNogs_egg6.tsv') as arc_levels:
    for line in arc_levels:
        if line.strip() and line.startswith('2157'):
            line = line.strip()
            info = line.split('\t')
            mems = info[3].split()

print (len(mems))

for m in ori_arc_sp:
    if m not in mems:
        print (m)
        lin = ncbi.get_lineage(m)
        print (lin)