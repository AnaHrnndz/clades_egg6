#!/bin/bash
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=main
#SBATCH -o /home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree/ete3_build_bact_sp_tree_%J.out
#SBATCH -e /home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree/ete3_build_bact_sp_tree_%J.err
#SBATCH -t 07-00:00

module unload ETE

#python2.7
ete3 build --cpu 40 -w clustalo_default-trimal01-none-none \
    -m cog_all-alg_concat_default-fasttree_default \
    -o /home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree \
    --clearall -a /home/plaza/projects/eggnog6/sp_tree/COG_files/Bact/Bact_COGs_seqs.faa \
    --cogs /home/plaza/projects/eggnog6/sp_tree/COG_files/Bact/COG_table.tsv \
    --noimg --spname-delimiter @

