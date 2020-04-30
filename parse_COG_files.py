from ete3 import SeqGroup
import sys
import os

f_ = sys.argv[1]
alg = SeqGroup(f_)
orginal_name_COG = os.path.basename(f_)

parse_name_COG = open('/home/plaza/projects/eggnog6/sp_tree/COG_files/Bact/clean_'+orginal_name_COG, 'w')
cog_tab = open('/home/plaza/projects/eggnog6/sp_tree/COG_files/Bact/COG_table.tsv', 'a')

sp_list_egg = []
with open('/home/plaza/projects/eggnog6/ncbi_tree/bact_sp_tree.tsv') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            sp = int(line.strip())
            sp_list_egg.append(sp)

eggnog62progenomes2 = {}
with open('/home/plaza/projects/eggnog6/raw_data/eggnogv6_2_progenomesv2.tsv') as tablefn:
    for line in tablefn:
        if line.strip() and not line.startswith("#"):
            info = line.split()
            taxid_egg6 = info[0].strip()
            taxid_pro2 = info[1].strip()
            eggnog62progenomes2[int(taxid_egg6)] = int(taxid_pro2)
            
sp_list = []
for sp in sp_list_egg:
    sp_pro = eggnog62progenomes2[sp]
    sp_list.append(int(sp_pro))

seq_list=[]
taxid_in_COG = []
for num, (name, seq, _) in enumerate(alg):
    name_data = name.split()
    new_name = name_data[0]
    try:
        taxid, sample, seq_name = new_name.split('.', 2)
        if int(taxid) in sp_list and int(taxid) not in taxid_in_COG:
            taxid_in_COG.append(int(taxid))
            final_name = str(taxid+'@'+sample+'-'+seq_name)
            seq_list.append(final_name)
            parse_name_COG.write('>%s\n%s\n' %(final_name, seq))
    except:
        print(new_name)

parse_name_COG.close()

cog_tab.write('\t'.join(seq_list)+'\n')
cog_tab.close()
 
