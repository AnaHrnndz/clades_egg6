sp_list = []
with open('/scratch/plaza/projects/eggnog6/guide_tree_ncbi/final_taxid_list') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            sp = line.strip()
            sp_list.append(sp)

print('len sp list: ', len(sp_list))


nogs_dir = {}
out_F = open('/scratch/plaza/projects/eggnog6/scripts/newXNogs_clean.txt', 'w')
with open('/scratch/plaza/projects/eggnog6/scripts/newXNogs.txt') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            info = line.split('\t')
            level = info[0]
            taxids_raw = list(info[4].split())
            taxids_clean = []
            for taxid in taxids_raw:
                if taxid in sp_list:
                    taxids_clean.append(taxid)
            nogs_dir[level] = taxids_clean
            if len(taxids_clean) != 0:
                out_F.write('\t'.join([level,info[1],info[2],info[3],(' '.join(nogs_dir[level])),'\n']))

out_F.close()


