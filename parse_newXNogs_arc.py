from ete3 import NCBITaxa

ncbi = NCBITaxa('/home/plaza/projects/eggnog6/ncbi_taxa_db/ete3_ncbitaxa_eggv5')

sp_list = []
with open('/home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/arc_final_taxid_list') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            sp = line.strip()
            sp_list.append(sp)

print('len sp list: ', len(sp_list))


nogs_dir = {}
out_F = open('/home/plaza/projects/eggnog6/clades/arc/newXNogs_arc.txt', 'w')
with open('/home/plaza/projects/eggnog6/raw_data/newXNogs.txt') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            info = line.split('\t')
            level = info[0]
            try:
                lineage = ncbi.get_lineage(level)
            except:
                name_class = info[2].split(':')[1]
                name2taxid = ncbi.get_name_translator([name_class])
                print(level, name_class,name2taxid)
                #level = name2taxid[name_class]
                #lineage = ncbi.get_lineage(level)
            if '2157' in lineage:
                taxids_raw = list(info[4].split())
                taxids_clean = []
                for taxid in taxids_raw:
                    if taxid in sp_list:
                        taxids_clean.append(taxid)
                nogs_dir[level] = taxids_clean
                if len(taxids_clean) != 0:
                    out_F.write('\t'.join([level,info[1],info[2],info[3],(' '.join(nogs_dir[level])),'\n']))

out_F.close()


