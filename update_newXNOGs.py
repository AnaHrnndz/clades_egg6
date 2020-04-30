from ete3 import NCBITaxa

ncbi = NCBITaxa('/home/plaza/projects/eggnog6/ncbi_taxa_db/taxa.sqlite')

eggv6_sp_list = []
with open('/home/plaza/projects/eggnog6/raw_data/eggNOGv6.species.list') as egg_list:
    for line in egg_list:
        if line.strip() and not line.startswith("#"):
            info = line.split('\t')
            sp = int(info[0])
            eggv6_sp_list.append(sp)
print (len(eggv6_sp_list))

#total_newXNogs_eggnogv6 = open('/home/plaza/projects/eggnog6/total_newXNogs_egg6.tsv', 'w')
#arc_newXNogs_v6 = open('/home/plaza/projects/eggnog6/arc_newXNogs_egg6.tsv', 'w')
bact_newXNogs_v6 = open('/home/plaza/projects/eggnog6/3_bact_newXNogs_egg6.tsv', 'w')
#euk_newXNogs_v6 = open('/home/plaza/projects/eggnog6/2_euk_newXNogs_egg6.tsv', 'w')
update_table = open('/home/plaza/projects/eggnog6/update_levels_eggnog6.tsv', 'w')


sp_parsed = list()
with open('/home/plaza/projects/eggnog6/raw_data/newXNogs.txt') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            info = line.split('\t')
            level = info[0]
            print (level)

            try:
                level = str(level)
                rank = ncbi.get_rank([level])
                tax2name = ncbi.get_taxid_translator([level])
                name = tax2name[int(level)]
                update_taxid = level
                lineage = ncbi.get_lineage(update_taxid)
                update_table.write(level+'\t'+name+'\t'+update_taxid+'\n')
    
            except:
                level = str(level)
                tax2name = ncbi.get_taxid_translator([level])
                name = tax2name[int(level)]
                name2newtaxid = ncbi.get_name_translator([name])
                update_taxid = name2newtaxid[str(name)][0]
                rank = ncbi.get_rank([update_taxid])
                lineage = ncbi.get_lineage(update_taxid)
                update_table.write(level+'\t'+str(name)+'\t'+str(update_taxid)+'\n')

            #descendants = ncbi.get_descendant_taxa(name)
            
            if not update_taxid in sp_parsed:

                member_taxid = []
                for tax in eggv6_sp_list:
                    lin_mem = ncbi.get_lineage(tax)
                    if int(update_taxid) in lin_mem:
                        member_taxid.append(tax)

                # if 2759 in lineage:
                    # euk_newXNogs_v6.write(str(update_taxid)+'\t'+rank[int(update_taxid)]+':'+name+'\t'+str(len(member_taxid))+'\t'+(' '.join(map(str, (member_taxid))))+'\n')
                

                if 2 in lineage:
                    bact_newXNogs_v6.write(str(update_taxid)+'\t'+rank[int(update_taxid)]+':'+name+'\t'+str(len(member_taxid))+'\t'+(' '.join(map(str, (member_taxid))))+'\n')
                
                sp_parsed.append(sp)
            # if 2157 in lineage:
            #     arc_newXNogs_v6.write(str(update_taxid)+'\t'+rank[int(update_taxid)]+':'+name+'\t'+str(len(member_taxid))+'\t'+(' '.join(map(str, (member_taxid))))+'\n')

            #print (str(update_taxid)+'\t'+rank[int(update_taxid)]+':'+name+'\t'+str(len(member_taxid))+'\t'+(' '.join(map(str, (member_taxid)))))
            #total_newXNogs_eggnogv6.write(str(update_taxid)+'\t'+rank[int(update_taxid)]+':'+name+'\t'+str(len(member_taxid))+'\t'+(' '.join(map(str, (member_taxid))))+'\n')


#euk_newXNogs_v6.close()
bact_newXNogs_v6.close()
# arc_newXNogs_v6.close()
# total_newXNogs_eggnogv6.close()