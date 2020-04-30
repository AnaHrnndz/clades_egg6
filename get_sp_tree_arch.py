##31-01-2020
##script created for building sp tree for freeze12 with ete3-NCBITaxa

from ete3 import NCBITaxa
ncbi = NCBITaxa('/home/plaza/projects/eggnog6/ncbi_taxa_db/taxa.sqlite')
#ncbi = NCBITaxa('/home/plaza/projects/eggnog6/ncbi_taxa_db/taxa.sqlite')

##Update ncbi taxonomy database with old version for progenomes2
#ncbi.update_taxonomy_database( taxdump_file='/home/plaza/projects/eggnog6/ref_tree_ncbi/taxdump/taxdump.tar.gz')


sci_name2taxid_ori = {}
sci_names_list = []
taxid_ori_list = []
with open('/home/plaza/projects/eggnog6/raw_data/Archaea_eggNOGv6.species.list') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            line_data = line.split('\t')
            taxid = line_data[0].strip()
            sci_name = line_data[1].strip()
            sci_name2taxid_ori[sci_name] = taxid
            sci_names_list.append(sci_name)
            taxid_ori_list.append(taxid)

print('len sciname2taxid ', len(sci_name2taxid_ori))
print('len sci names list', len(sci_names_list))


taxid_ncbi = []
for name in sci_names_list:
    name2taxid_ncbitaxa = ncbi.get_name_translator([name])
    if len(name2taxid_ncbitaxa) == 0:
        print(name)
        taxid = sci_name2taxid_ori[name]
        taxid_ncbi.append(taxid)


name2taxid_ncbitaxa = ncbi.get_name_translator(sci_names_list)
print('len name2taxid ',len(name2taxid_ncbitaxa))


for name, taxid in name2taxid_ncbitaxa.items():
    taxid_ncbi.append(taxid[0])

print ('len taxid ncbi', len(taxid_ncbi))
tree = ncbi.get_topology(taxid_ncbi)
print('len tree ', len(tree))

leaf_name_list = []
for leaf in tree:
    taxid_ncbi = leaf.name
    taxid2sciname = ncbi.get_taxid_translator([taxid_ncbi])

    try:
        sci_name_ncbi = taxid2sciname[int(taxid_ncbi)]
        leaf.name = sci_name2taxid_ori[sci_name_ncbi]
    except:
        leaf.name = taxid_ncbi
    
    leaf_name_list.append(leaf.name)

print('len tree ', len(tree))

out_f = open('/home/plaza/projects/eggnog6/ncbi_tree/arch_sp_tree.tsv', 'w')
for leaf in leaf_name_list:
    out_f.write(leaf+'\n')
out_f.close()

miss_names = []
for name in taxid_ori_list:
    if name in leaf_name_list:
        continue
    else:
        miss_names.append(name)
print(len(miss_names), miss_names)



tree.write(format = 1 , outfile = '/home/plaza/projects/eggnog6/ncbi_tree/Archaea_ncbi.nw')
        #print('miss sp: ', name)


# for leaf in tree:
#     leaf_sci_name = leaf.name
#     ori_sci_name = sci_name2taxid_ori[leaf_sci_name]
#     leaf.name = ori_sci_name


# miss_names_2 = []
# for name in sci_names_list:
#     if name in leaf_name_list:
#         continue
#     else:
#         miss_names_2.append(name)
# print(len(miss_names_2))

# sp_list = []
# with open('/home/plaza/projects/eggnog6/ref_taxid_list.txt') as tablefn:
#     for line in tablefn:
#         if line.strip() and not line.startswith("#"):
#             sp = int(line.split('>')[1])
#             sp_list.append(sp)

# print('len sp list: ', len(sp_list))

# tree = ncbi.get_topology(sp_list)
# print ('num leafs:', len(tree))

# n_names = set()
# for n in tree.traverse():
#     n_names.add(str(n.name))
# print('number nodes: ', len(n_names))

# miss_sp = []
# miss_sp_t = {}
# for sp in sp_list:
#     if str(sp) not in n_names:
#         miss_sp.append(str(sp))
#         new_name = good2new[str(sp)]
#         miss_sp_t[new_name] = str(sp)
# print('miss sp: ', len(miss_sp))
# miss_1=open('/home/plaza/projects/eggnog6/guide_tree_ncbi/miss_1.tsv','w')
# miss_1.write('\n'.join(miss_sp))
# miss_1.close()


# n_names_2 = set()
# for n in tree.traverse():
#     if n.name == str(1458425):
#         n.name = str(1458425)
#         n_names_2.add(str(n.name))
#     elif n.name == str(656366):
#         n.name = str(656366)
#         n_names_2.add(str(n.name))
#     elif n.name in miss_sp_t:
#         new_name = str(n.name)
#         good_name = miss_sp_t[new_name]
#         n.name = str(good_name)
#         n_names_2.add(str(n.name))
#     else:
#         n_names_2.add(str(n.name))

# miss_sp_2 =  []
# for sp in sp_list:
#     if str(sp) not in n_names_2:
#         miss_sp_2.append(str(sp))
# print('miss sp 2: ', len(miss_sp_2))
# miss_2=open('/home/plaza/projects/eggnog6/guide_tree_ncbi/miss_2.tsv','w')
# miss_2.write('\n'.join(miss_sp_2))
# miss_2.close()


# tree.write(format=1, outfile="/home/plaza/projects/eggnog6/guide_tree_ncbi/freeze12_representative_sp_tree_oldDB.nw")